import os
import sys
import h5py
import pickle
import shutil
import numpy as np
import subprocess
from fluent_class.fluent_simulation import fluentSimulation

from ofClass.fluned_tool_launchers import launch_volume_func_object
from ofClass.fluned_tool_launchers import launch_grad_func_object
from .util import mod


from ofClass.ofClass import SimulationOF
from ofClass.fluned_vtk_utils import write_cartesian_vtk


class flunedCase:
    """
    FLUNED simulation class
    """

    def __init__(self, fluned_path=""):
        """
        initialize flunedCase
        """
        self.decay_constant = 0
        self.activation_normalization = 0
        self.activation_dataset = ""
        self.activation_dataset_error = ""
        self.activation_const = 0
        self.activation_file = ""
        self.fv_scheme = "stable"
        self.fluned_path = fluned_path
        self.case = os.path.split(fluned_path)[1]
        self.num_internal_cells = 0

    def generate_case(self, arg_dict):
        """
        fill simulation parameters
        """

        self.case = arg_dict["case"]
        # set decay constant

        # default value
        if "simulation_type" not in arg_dict:
            arg_dict["simulation_type"] = "single-isotope"

        if arg_dict["simulation_type"] == "single-isotope":
            if "isotope" in arg_dict:
                isotope = arg_dict["isotope"].lower().replace("-", "")
                if isotope not in ["n16", "o19", "n17", "f20", "custom"]:
                    raise ValueError("isotope not recognized")
                self.isotope = isotope
            else:
                self.isotope = "custom"

            if "decay_constant" in arg_dict:
                self.decay_constant = float(arg_dict["decay_constant"])

            if "activation_dataset" in arg_dict:
                self.activation_dataset = arg_dict["activation_dataset"]

            if "activation_dataset_error" in arg_dict:
                self.activation_dataset_error = arg_dict["activation_dataset_error"]

            if "activation_file" in arg_dict:
                self.activation_file = os.path.normcase(arg_dict["activation_file"])

        if arg_dict["simulation_type"] == "openmc-multi":
            import openmc.data

            self.isotope = arg_dict["isotope"]
            isotope_string = self.isotope[0].upper() + self.isotope[1:].lower()

            self.decay_constant = openmc.data.decay_constant(self.isotope)

            # mapping of the reaction rates to the vtk file

            reaction_rate_file_path = arg_dict["activation_file"]
            reaction_rate_file_dir = os.path.dirname(arg_dict["activation_file"])

            activation_file_name = f"reaction_rate_{self.isotope.upper()}.vtk"
            activation_file_path = os.path.join(
                reaction_rate_file_dir, activation_file_name
            )
            activation_dataset = "reaction_rate_m3"

            self.activation_dataset = activation_dataset

            with h5py.File(reaction_rate_file_path, "r") as f:
                neutron_fluxes = pickle.loads(f["fluxes"][...].tobytes())
                openmc_strenght = pickle.loads(f["openmc_strength"][...].tobytes())
                isotope_densities = pickle.loads(f["isotope_density"][...].tobytes())
                mesh_widths_cm = pickle.loads(f["mesh_width"][...].tobytes())
                mesh_dimensions = pickle.loads(f["mesh_dimension"][...].tobytes())
                mesh_lower_left_cm = pickle.loads(f["mesh_lower_left"][...].tobytes())
                micro_xs_data = f["data"][:]

            # adapt to meters and pyvista format
            mesh_lower_left_m = [coord * 0.01 for coord in mesh_lower_left_cm]
            mesh_widths_m = [width * 0.01 for width in mesh_widths_cm]
            mesh_dimensions_pv = [dim + 1 for dim in mesh_dimensions]

            macro_xs_channels = []

            for channel in arg_dict["openmc_depletion_channels"]:
                idxs = channel["channel_index"]

                parent = channel["parent_nuclide"]

                macro_xs_channel = (
                    micro_xs_data[:, idxs[0], idxs[1]] * isotope_densities[parent]
                )

                macro_xs_channels.append(macro_xs_channel)

            total_rr_m3 = (
                np.multiply(np.sum(macro_xs_channels, axis=0), neutron_fluxes)
            ) * openmc_strenght
            # NB no scaling factor is hardcoded to go from cm-3 to m-3.
            # either modify the input data or the fluned input file in the appropriate option

            # print("debug - total rr :", sum(total_rr_m3))

            write_cartesian_vtk(
                activation_file_name,
                mesh_dimensions_pv,
                mesh_widths_m,
                mesh_lower_left_m,
                total_rr_m3,
                activation_dataset,
            )

            self.activation_file = activation_file_path

        if "activation_const" in arg_dict:
            self.activation_const = float(arg_dict["activation_const"])

        if "activation_normalization" in arg_dict:
            self.activation_normalization = float(arg_dict["activation_normalization"])

        if "fv_scheme" in arg_dict:
            if arg_dict["fv_scheme"].lower() in ["stable", "accurate"]:
                self.fv_scheme = arg_dict["fv_scheme"]
            else:
                raise ValueError("fv_scheme must be 'stable' or 'accurate'")

        if "cfd_type" not in arg_dict:
            raise ValueError("ERROR: cfd type not specified in input file")
        else:
            self.cfd_type = arg_dict["cfd_type"].lower()

        if "time_treatment" not in arg_dict:
            raise ValueError("ERROR: type of time treatment not specified")
        elif arg_dict["time_treatment"].lower() not in ["steadystate", "transient"]:
            err_string = "ERROR: time treatment must be 'steadyState' or 'transient'"
            raise ValueError(err_string)
        else:
            self.time_treatment = arg_dict["time_treatment"].lower()

        if self.cfd_type in ["fluent-h5-multi", "fluent-multi"]:
            if "fluent_fluid_region_name" not in arg_dict:
                print("ERROR: name of the fluid region to extract not")
                print("specified! use parameter FLUENT_FLUID_REGION_NAME")
                raise ValueError()
            else:
                self.fluent_fluid_region_name = arg_dict["fluent_fluid_region_name"]

        self.molecular_diffusion = float(arg_dict["molecular_diffusion"])
        self.schmidt_number = float(arg_dict["schmidt_number"])
        self.inlet_conc = float(arg_dict["inlet_conc"])

        self.cfd_path = os.path.normcase(arg_dict["cfd_path"])

        if not os.path.isdir(self.cfd_path):
            raise OSError(f"Folder not found: {self.cfd_path}")

        self.fluned_path = os.path.join(self.cfd_path, self.case)

        # try:
        #     Path(self.fluned_path).mkdir(parents=True)
        # except FileExistsError:
        #     pass

        return

    def generate_fluned_files(self):
        """
        this function initialize the cfd class, after the initial openfoam or
        fluent simulation has been processed
        """

        self.cfd_simulation.post_process_openfoam_simulation()

        self.fluned_simulation = SimulationOF(self.fluned_path)

        self.fluned_simulation.create_case_folder()

        self.copy_last_phi()
        self.copy_last_u()
        self.copy_last_nut()
        self.copy_poly_mesh()

        self.fluned_simulation.patches = self.cfd_simulation.patches
        self.fluned_simulation.n_internal_cells = self.cfd_simulation.n_internal_cells

        self.fluned_simulation.generate_system_files(
            self.time_treatment, self.fv_scheme
        )

        self.fluned_simulation.generate_constant_file(
            self.isotope,
            self.molecular_diffusion,
            self.decay_constant,
            self.schmidt_number,
        )

        self.fluned_simulation.generate_zero_t(self.inlet_conc)
        self.fluned_simulation.generate_zero_ta()
        self.fluned_simulation.generate_zero_td(self.inlet_conc)
        self.fluned_simulation.generate_zero_tr(self.time_treatment)

        launch_volume_func_object(self.fluned_path)

        self.fluned_simulation.assign_activation_rates(
            self.activation_file,
            self.activation_const,
            self.activation_dataset,
            self.activation_dataset_error,
            self.activation_normalization,
        )
        self.fluned_simulation.generate_source_file()
        self.fluned_simulation.generate_tr_source_file(self.time_treatment)

        return

    def parse_fluned_simulation(self):
        """
        this function gets the data of a completed FLUNED simulation
        it initialize a SimulationOF object and parse its data
        """

        self.fluned_simulation = SimulationOF(self.fluned_path)

        self.fluned_simulation.post_process_openfoam_simulation()
        self.fluned_simulation.post_process_fluned_simulation()

        return

    def initialize_cfd_class(self):
        """
        this function initialize the cfd class, after the initial openfoam or
        fluent simulation has been processed
        """

        self.cfd_simulation = SimulationOF(self.cfd_path)
        # self.cfd_simulation.create_case_folder()

        return

    def create_case_folders_fluent(self):
        """
        This function creates additional folders required for the FLUENT simulation converted
        into the OpenFOAM format
        """

        case_folder = self.cfd_path

        # T=1 folder
        one_folder = os.path.join(case_folder, "1")
        if not os.path.exists(one_folder):
            os.mkdir(one_folder)

        return

    def convert_fluent_to_openfoam(self):
        """
        it writes the fluent data in the openfoam format
        """

        self.fluent_simulation.convert_to_openfoam(self.cfd_path)

        return

    def parse_fluent_simulation(self):
        """
        this function initialize the fluent class
        """

        self.fluent_simulation = fluentSimulation(
            self.cfd_path, self.fluent_fluid_region_name
        )

        self.fluent_simulation.parse_multidataset_h5_file()

        return

    def copy_last_phi(self):
        """
        this function look for the last phi file in the cfd folder
        if in the cfd simulation the flow is in m3/s it just copies it to the
        fluned folders. If the flows is in kg/s it converts it to m3/s
        """
        # print("debug: copying phi file")
        last_time = self.cfd_simulation.last_time
        target_file = os.path.join(self.fluned_path, "0", "phi")
        origin_file = os.path.join(self.cfd_path, str(last_time), "phi")

        if self.cfd_simulation.volumetric_flag:
            shutil.copyfile(origin_file, target_file)
            # print("debug: copying phi file action")

        else:
            print("converting phi values to volumetric flow ...")
            self.cfd_simulation.convert_phi_to_volumetric(target_file)
            # print("debug: converting to volumetric")

        return

    def copy_last_u(self):
        """
        this function look for the last U file in the cfd folder
        """

        last_time = self.cfd_simulation.last_time
        target_file = os.path.join(self.fluned_path, "0", "U")
        origin_file = os.path.join(self.cfd_path, str(last_time), "U")
        shutil.copyfile(origin_file, target_file)

        return None

    def copy_last_nut(self):
        """this function look for the last nut file in the cfd folder,
        if it does not exists it means the simulation is laminar"""

        last_time = self.cfd_simulation.last_time

        targetFile = os.path.join(self.fluned_path, "0", "nut")
        originFile = os.path.join(self.cfd_path, str(last_time), "nut")
        checkFile = os.path.isfile(originFile)
        if checkFile:
            shutil.copyfile(originFile, targetFile)

        return

    def copy_poly_mesh(self):
        """this function copy the poly Mesh files in the FLUNED folder"""

        sourceFolder = os.path.join(self.cfd_path, "constant", "polyMesh")
        targetFolder = os.path.join(self.fluned_path, "constant", "polyMesh")

        # shutil do not allow to copy into an existing folder
        dirCheck = os.path.isdir(targetFolder)
        if dirCheck:
            shutil.rmtree(targetFolder)

        shutil.copytree(sourceFolder, targetFolder)

        return

    def launch_solver(self):
        """
        this function launches the scalar calculation
        """

        print("launching FLUNED solver ...")
        current_folder = os.getcwd()
        os.chdir(self.fluned_path)
        launch_calc_string = "FLUNED-solver".split()
        with open("simulation_log", "a", encoding="utf-8") as outfile:
            subprocess.Popen(launch_calc_string, stdout=outfile).wait()
        os.chdir(current_folder)

        return

    def launch_check_mode(self):
        """
        this function launches some functional objects and read its results for some additional checks
        """

        launch_grad_func_object(self.fluned_path)
        self.fluned_simulation.read_velocities()

        return

    def generate_cartesian_radiation_source_model(self):
        """
        this function sample the fluned results into a cartesian mesh
        """

        self.fluned_simulation.calculate_cartesian_sampling_coordinates(
            self.source_sampling_resolution_cm
        )
        self.fluned_simulation.sample_source_to_cartesian_mesh(
            self.source_sampling_dataset
        )

        self.fluned_simulation.write_sampled_cartesian_source_vtk()

        return

    def write_results(self, arguments):
        """
        this function generate the final results
        """

        if self.fluned_simulation.time_treatment == "steadystate":
            self.write_summary_steady(arguments)
        else:
            raise NotImplementedError("ERROR transient summary not implemented yet")
            # self.write_summary_transient(arguments)

        return None

    def write_summary_transient(self, arguments):
        """
        this function is not ready yet after the large 2025 refactoring
        """
        # sumFile = os.path.join(self.results_folder, "SUMMARY.csv")
        #
        # inletAtoms = 0
        # inletFlow = 0
        # outletAtoms = 0
        # outletFlow = 0
        #
        # patches = sorted(self.patches.values(), key=lambda x: x.face_id)
        #
        # for face in patches:
        #     if face["typeFile"] == "inlet":
        #         inletAtoms += abs(face["TFlowFileLast"])
        #         inletFlow += abs(face["phiFlowFileLast"])
        #
        #     if face["typeFile"] == "outlet":
        #         outletAtoms += abs(face["TFlowFileLast"])
        #         outletFlow += abs(face["phiFlowFileLast"])
        #
        # totalInletConc = inletAtoms / inletFlow
        # totalInActivity = totalInletConc * self.decay_constant
        #
        # totalOutletConc = outletAtoms / outletFlow
        # totalOutActivity = totalOutletConc * self.decay_constant
        #
        # cellActivity = np.zeros(len(self.t_scalar))
        #
        # for i in range(len(self.t_scalar)):
        #     if self.t_scalar[i] > 0:
        #         cellActivity[i] = (
        #             self.TScalar[i] * self.Volumes[i] * self.decay_constant
        #         )
        #     else:
        #         cellActivity[i] = 0
        #
        # totActivity = sum(cellActivity)
        # avgActivity = totActivity / sum(self.volumes)
        #
        # averageVolume = sum(self.volumes) / len(self.volumes)
        #
        # negCheck = any(n < 0 for n in self.t_scalar)
        #
        # with open(sumFile, "w") as fw:
        #     fw.write("FLUNED SIMULATION SUMMARY\n")
        #     fw.write("CASE,{},\n".format(self.case))
        #     fw.write("TRANSIENT SIMULATION\n")
        #     fw.write("N ELEMENTS,{},\n".format(len(self.t_scalar)))
        #     fw.write("ISOTOPE,{},\n".format(self.isotope.upper()))
        #     fw.write("DECAY CONSTANT,{:e},\n".format(self.decay_constant))
        #     fw.write("MOL DIFFUSION,{:e},\n".format(self.molecular_dif))
        #     fw.write("TURB SCHMIDT N,{:f},\n".format(self.schmidt_number))
        #     fw.write("\n")
        #     fw.write("\n")
        #     fw.write("QUALITY\n")
        #     fw.write("AVG VOL [m3],{:e},\n".format(averageVolume))
        #
        #     if negCheck:
        #         fw.write("WARNING some elements are negative")
        #
        #     if arguments.cdgs:
        #         fw.write("\n")
        #         fw.write("\n")
        #         fw.write("SOURCE SAMPLING\n")
        #         fw.write("SAMPLING RESOLUTION [m],{:f},\n".format(self.precision))
        #         fw.write(
        #             "SAMPLED VOXELS [#],{:d},\n".format(
        #                 self.xInts * self.yInts * self.zInts
        #             )
        #         )
        #         fw.write(
        #             "VTK EMISSION RATE [#/s],{:e},\n".format(self.originalEmissionRate)
        #         )
        #         sampString = "SAMPLED EMISSION RATE (UNSCALED) [#/s],{:e},\n"
        #         fw.write(sampString.format(self.unscaledEmissionRate))
        #
        #     fw.write("\n")
        #     fw.write("\n")
        #     fw.write("ACTIVATION\n")
        #     fw.write("INLET ATOMS FINAL [#/s],{:.5e},\n".format(inletAtoms))
        #     fw.write("OUTLET ATOMS FINAL [#/s],{:.5e},\n".format(outletAtoms))
        #     fw.write("TOT IN ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(totalInActivity))
        #     fw.write(
        #         "TOT OUT ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(totalOutActivity)
        #     )
        #     fw.write("TOT CASE ACTIVITY FINAL [Bq],{:.5e},\n".format(totActivity))
        #     fw.write("TOT AVG ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(avgActivity))
        #     if totalInActivity != 0:
        #         normAvgActivity = avgActivity / totalInActivity
        #         stri = "INLET-NORMALIZED AVG ACTIVITY FINAL [Bq/m3],{:.5e},\n"
        #         fw.write(stri.format(normAvgActivity))
        #
        #     if inletAtoms != 0:
        #         reductionRate = outletAtoms / inletAtoms
        #         fw.write("OUT/IN RATIO FINAL,{:.5e},\n".format(reductionRate))
        #
        #     fw.write("\n")
        #     fw.write("\n")
        #     fw.write("FACES\n")
        #
        #     for face in patches:
        #         if face.face_type != "wall":
        #             fw.write(face.face_id)
        #             fw.write("\n")
        #             fw.write("TYPE,{},\n".format(face.face_type))
        #             fw.write("AREA [m2],{:.5e},\n".format(face["areaFile"]))
        #             fw.write(
        #                 "FLUID FLOW FINAL [m3/s],{:.5e},\n".format(
        #                     face["phiFlowFileLast"]
        #                 )
        #             )
        #             fw.write(
        #                 "ATOM FLOW FINAL [#/s],{:.5e},\n".format(face["TFlowFileLast"])
        #             )
        #             fw.write(
        #                 "ATOM CONC FINAL [#/m3],{:.5e},\n".format(face["avTFileLast"])
        #             )
        #             fw.write(
        #                 "SPECIFIC ACTIVITY FINAL [Bq/m3],{:.5e},\n".format(
        #                     abs(face["avTFileLast"]) * self.decay_constant
        #                 )
        #             )
        #             if face["typeFile"] == "outlet":
        #                 fluxFrac = abs(face["phiFlowFileLast"] / outletFlow)
        #                 fw.write(
        #                     "AVG RTD RES T [s],{:.5f},\n".format(
        #                         abs(face["rtdResTime"])
        #                     )
        #                 )
        #                 fw.write(
        #                     "RTD FACE RED RATE,{:.5f},\n".format(
        #                         fluxFrac * (face["rtdDecRate"])
        #                     )
        #                 )
        #             fw.write("\n")
        #
        # for face in self.facesPost:
        #     if face["typeFile"] == "wall":
        #         continue
        #
        #     # write the flowing atoms at inlet-outlet
        #     faceSummary = (
        #         "face-atom-flow-" + face["faceID"] + "-" + face["typeFile"] + ".csv"
        #     )
        #     sumFile1 = os.path.join(self.results_folder, faceSummary)
        #     with open(sumFile1, "w") as fw:
        #         for t, c in zip(face["timeFile"], face["TFlowFile"]):
        #             fw.write("{:.3f},{:.5e},\n".format(t, c))
        #
        #     # write the concentration at inlet-outlet
        #     faceSummary = (
        #         "face-conc-" + face["faceID"] + "-" + face["typeFile"] + ".csv"
        #     )
        #     sumFile1 = os.path.join(self.results_folder, faceSummary)
        #     with open(sumFile1, "w") as fw:
        #         for t, c in zip(face["timeFile"], face["avTFile"]):
        #             fw.write("{:.3f},{:.5e},\n".format(t, c))
        #
        #     # write the specific activity at inlet-outlet
        #     faceSummary = (
        #         "face-specific-activity-"
        #         + face["faceID"]
        #         + "-"
        #         + face["typeFile"]
        #         + ".csv"
        #     )
        #     sumFile1 = os.path.join(self.results_folder, faceSummary)
        #     with open(sumFile1, "w") as fw:
        #         for t, c in zip(face["timeFile"], face["avTFile"]):
        #             fw.write("{:.3f},{:.5e},\n".format(t, c * self.decay_constant))
        #
        #     # write the specific time-conc at inlet-outlet
        #     faceSummary = (
        #         "face-fictitious-time-"
        #         + face["faceID"]
        #         + "-"
        #         + face["typeFile"]
        #         + ".csv"
        #     )
        #     sumFile1 = os.path.join(self.results_folder, faceSummary)
        #     with open(sumFile1, "w") as fw:
        #         for t, c in zip(face["timeFile"], face["avTrFile"]):
        #             fw.write("{:.3f},{:.5e},\n".format(t, c))
        #
        #     # write the RTD at inlet-outlet
        #     faceSummary = (
        #         "face-RTD-raw-" + face["faceID"] + "-" + face["typeFile"] + ".csv"
        #     )
        #     sumFile1 = os.path.join(self.results_folder, faceSummary)
        #     with open(sumFile1, "w") as fw:
        #         for t, g in zip(face["timeFile"], face["avTrGrad"]):
        #             fw.write("{:.3f},{:.5e},\n".format(t, g))
        #
        # return
        #

    def write_summary_steady(self, arguments):
        """
        This function writes the summary file in the RESULTS/ folder
        """

        summary_file = os.path.join(
            self.fluned_simulation.results_folder, "SUMMARY.csv"
        )

        results_sim = self.fluned_simulation

        inlet_atoms = abs(results_sim.total_inlet_t_atoms)
        inlet_activity = abs(results_sim.inlet_td_conc_atoms_m3[-1])

        outlet_activity = results_sim.outlet_t_conc_atoms_m3[-1]
        outlet_atoms = results_sim.total_outlet_t_atoms

        tot_activity = results_sim.total_isotope_activity
        avg_activity = results_sim.total_average_isotope_concentration

        averageVolume = sum(results_sim.volumes) / len(results_sim.volumes)

        faces_post = sorted(results_sim.patches.values(), key=lambda x: x.face_id)

        if arguments.check:
            lbda = self.decay_constant
            ln2 = np.log(2)

            totS = sum(results_sim.t_scalar)

            # parameter 2 based on transit times

            transitTimes = [
                (vol ** (1 / 3)) / mod(v)
                for vol, v in zip(results_sim.volumes, results_sim.velocities)
            ]

            averageTransitTime = sum(transitTimes) / len(transitTimes)

            qualityPar2 = [
                t * c * lbda / (ln2 * totS)
                for t, c in zip(transitTimes, results_sim.t_scalar)
            ]

            avgMeshQualParameter2 = sum(qualityPar2) / len(qualityPar2)

            # parameter 3 based on scalar gradient

            gradients = [mod(g) for g in results_sim.gradients]

            gradientDist = [
                g * (vol ** (1 / 3)) for vol, g in zip(results_sim.volumes, gradients)
            ]

            avgGradDist = sum(gradientDist) / len(gradientDist)

            qualityPar3 = [
                gd * c / totS for gd, c in zip(gradientDist, results_sim.t_scalar)
            ]

            avgMeshQualParameter3 = sum(qualityPar3) / len(qualityPar3)

        negCheck = any(n < 0 for n in results_sim.t_scalar)

        with open(summary_file, "w") as fw:
            fw.write("FLUNED SIMULATION SUMMARY\n")
            fw.write("CASE,{},\n".format(results_sim.case))
            fw.write("STEADY STATE SIMULATION\n")
            fw.write("N ELEMENTS,{},\n".format(results_sim.n_internal_cells))
            fw.write("ISOTOPE,{},\n".format(results_sim.isotope.upper()))
            fw.write("DECAY CONSTANT,{:e},\n".format(results_sim.decay_constant))
            fw.write("MOL DIFFUSION,{:e},\n".format(results_sim.molecular_diffusion))
            fw.write("TURB SCHMIDT N,{:f},\n".format(results_sim.schmidt_number))
            fw.write("\n")
            fw.write("\n")
            fw.write("QUALITY\n")
            fw.write("AVG VOL [m3],{:e},\n".format(averageVolume))
            if arguments.check:
                fw.write("AVG CELL TRNST T [s],{:f},\n".format(averageTransitTime))
                fw.write("AVG CELL Q2,{:e},\n".format(avgMeshQualParameter2))
                fw.write("AVG CELL GRADxLEN,{:e},\n".format(avgGradDist))
                fw.write("AVG CELL Q3,{:e},\n".format(avgMeshQualParameter3))
            if negCheck:
                fw.write("WARNING some elements are negative")
            if arguments.cdgs:
                fw.write("\n")
                fw.write("\n")
                fw.write("SOURCE SAMPLING\n")
                fw.write(
                    "SAMPLING RESOLUTION [m],{:f},\n".format(
                        self.source_sampling_resolution_cm / 100
                    )
                )
                fw.write(
                    "SAMPLED VOXELS [#],{:d},\n".format(
                        results_sim.x_ints * results_sim.y_ints * results_sim.z_ints
                    )
                )
                fw.write(
                    "VTK EMISSION RATE [#/s],{:e},\n".format(
                        results_sim.total_isotope_emission_rate
                    )
                )
                sampString = "SAMPLED EMISSION RATE (UNSCALED) [#/s],{:e},\n"
                fw.write(sampString.format(results_sim.raw_sampled_total_emission_rate))
            fw.write("\n")
            fw.write("\n")
            fw.write("ACTIVATION\n")
            fw.write("INLET ATOMS [#/s],{:.5e},\n".format(inlet_atoms))
            fw.write("OUTLET ATOMS [#/s],{:.5e},\n".format(outlet_atoms))
            fw.write("TOT IN ACTIVITY [Bq/m3],{:.5e},\n".format(inlet_activity))
            fw.write("TOT OUT ACTIVITY[Bq/m3],{:.5e},\n".format(outlet_activity))
            fw.write("TOT CASE ACTIVITY[Bq],{:.5e},\n".format(tot_activity))
            fw.write("TOT AVG ACTIVITY[Bq/m3],{:.5e},\n".format(avg_activity))
            if inlet_activity != 0:
                norm_avg_activity = abs(avg_activity / inlet_activity)
                string = "INLET-NORMALIZED AVG ACTIVITY[Bq/m3],{:.5e},\n"
                fw.write(string.format(norm_avg_activity))

            if inlet_atoms != 0:
                reduction_rate = abs(outlet_atoms / inlet_atoms)
                fw.write("OUT/IN RATIO,{:.5e},\n".format(reduction_rate))

            fw.write("\n")
            fw.write("\n")
            fw.write("FACES\n")

            for face in faces_post:
                if face.face_type != "wall":
                    fw.write(face.face_id)
                    fw.write("\n")
                    fw.write("TYPE,{},\n".format(face.face_type))
                    fw.write("AREA [m2],{:.5e},\n".format(face.area_m2))
                    fw.write(
                        "FLUID FLOW [m3/s],{:.5e},\n".format(face.post_process_flow[-1])
                    )
                    fw.write(
                        "ATOM FLOW [#/s],{:.5e},\n".format(face.post_process_t_flow[-1])
                    )
                    fw.write(
                        "ATOM CONC [#/m3],{:.5e},\n".format(face.t_conc_atoms_m3[-1])
                    )
                    fw.write(
                        "SPECIFIC ACTIVITY [Bq/m3],{:.5e},\n".format(
                            abs(face.t_conc_atoms_m3[-1]) * results_sim.decay_constant
                        )
                    )
                    fw.write("AVG RES T [s],{:.5f},\n".format(abs(face.tr_conc[-1])))
                    fw.write("\n")

        return
