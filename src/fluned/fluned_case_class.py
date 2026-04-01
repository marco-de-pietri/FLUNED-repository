import pickle
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import numpy as np

from fluned.fluent_class.fluent_simulation import fluentSimulation
from fluned.isotopes.isotopes import load_isotopes
from fluned.ofoam_class.fluned_mesh_utils import (
    apply_roto_translation_to_vtk_grid,
    write_cartesian_vtk,
)
from fluned.ofoam_class.fluned_tool_launchers import (
    launch_grad_func_object,
    launch_volume_func_object,
)
from fluned.ofoam_class.ofoam_base import oFoamBase
from fluned.ofoam_class.ofoam_fluned import oFoamFluned

from .utils import mod


class flunedCase:
    """
    FLUNED simulation class - two class of jobs:
    - orchestrate the generation of the FLUNED case files from the initial CFD simulation results
    - orchestrate the post processing
    """

    def __init__(self, fluned_path: str | Path = ""):
        """
        initialize flunedCase
        """
        self.fluned_path = fluned_path

    def generate_case(self, arg_dict):
        """
        fill simulation parameters
        """

        self.case = str(arg_dict["case"])
        self.source_sampling_dataset = str(arg_dict["source_sampling_dataset"])
        self.source_sampling_resolution_cm = float(
            arg_dict["source_sampling_resolution_cm"]
        )

        if arg_dict["simulation_type"] == "single-isotope":
            # set isotope and decay constant in single-isotope simulations
            if "isotope" in arg_dict:
                isotope = arg_dict["isotope"].lower().replace("-", "")
                if isotope not in ["n16", "o19", "n17", "f20", "custom"]:
                    raise ValueError(
                        "isotope currently not allowed - limit to n16, n17, o19, f20"
                        "for more isotopes use openmc multi simulation"
                    )
                self.isotope = isotope

            if "decay_constant" in arg_dict and arg_dict["decay_constant"] != 0.0:
                self.decay_constant = float(arg_dict["decay_constant"])
            elif self.isotope != "custom":
                isotope_database = load_isotopes()
                iso_key = self.isotope  # already lowercase in this branch
                if iso_key not in isotope_database:
                    raise ValueError("ERROR isotope not found in the database")
                isotope_data = isotope_database[iso_key]
                self.decay_constant = isotope_data.decay_constant
            else:
                raise ValueError(
                    "neither the isotope nor the decay constant is specified, cannot \
                    proceed"
                )

            self.activation_dataset = arg_dict["activation_dataset"]
            self.activation_dataset_error = arg_dict["activation_dataset_error"]

            if "activation_file" in arg_dict:
                self.activation_file = arg_dict["activation_file"]

        if arg_dict["simulation_type"] == "openmc-multi":
            import openmc.data  # type: ignore[import]

            self.isotope = arg_dict["isotope"]

            self.decay_constant = openmc.data.decay_constant(self.isotope)

            # mapping of the reaction rates to the vtk file

            reaction_rate_file_path = Path(arg_dict["activation_file"])
            reaction_rate_file_dir = reaction_rate_file_path.parent

            activation_file_name = f"reaction_rate_{self.isotope.upper()}.vtk"
            activation_file_path = reaction_rate_file_dir / activation_file_name
            activation_dataset = "reaction_rate_m3"

            self.activation_dataset = activation_dataset
            self.activation_dataset_error = arg_dict["activation_dataset_error"]

            with h5py.File(reaction_rate_file_path, "r") as f:
                neutron_fluxes = pickle.loads(f["fluxes"][...].tobytes())
                isotope_densities = pickle.loads(f["isotope_density"][...].tobytes())
                mesh_widths_cm = pickle.loads(f["mesh_width"][...].tobytes())
                mesh_dimensions = pickle.loads(f["mesh_dimension"][...].tobytes())
                mesh_lower_left_cm = pickle.loads(f["mesh_lower_left"][...].tobytes())
                micro_xs_data = f["data"][:]

            # adapt to meters and pyvista format
            mesh_lower_left_m = [coord * 0.01 for coord in mesh_lower_left_cm]
            mesh_widths_m = [width * 0.01 for width in mesh_widths_cm]
            mesh_dimensions_pv = [dim + 1 for dim in mesh_dimensions]
            mesh_voxel_volume_cm3 = (
                mesh_widths_cm[0] * mesh_widths_cm[1] * mesh_widths_cm[2]
            )

            macro_xs_channels = []

            for channel in arg_dict["openmc_depletion_channels"]:
                idxs = channel["channel_index"]

                parent = channel["parent_nuclide"]

                macro_xs_channel = (
                    micro_xs_data[:, idxs[0], idxs[1]] * isotope_densities[parent]
                )

                macro_xs_channels.append(macro_xs_channel)

                # print(
                #     "debug: macro xs channel for {}: {}".format(
                #         channel["parent_nuclide"], macro_xs_channel
                #     )
                # )
                # print("parent nuclide: {}".format(parent))
                # print("mesh widths: {}".format(mesh_widths_m))
                # print("mesh dimensions: {}".format(mesh_dimensions_pv))
                # print("isotope densities: {}".format(isotope_densities[parent]))
                # print(
                #     "average microxs for channel {}: {}".format(
                #         channel["parent_nuclide"],
                #         np.mean(micro_xs_data[:, idxs[0], idxs[1]]),
                #     )
                # )
                # print("average flux: {}".format(np.mean(neutron_fluxes)))
                # print()

            # WARNING scaling factor is hardcoded to go from cm-3 to m-3.
            total_rr_m3 = (
                (np.multiply(np.sum(macro_xs_channels, axis=0), neutron_fluxes))
                / mesh_voxel_volume_cm3
            ) * 1e6

            write_cartesian_vtk(
                activation_file_name,
                mesh_dimensions_pv,
                mesh_widths_m,
                mesh_lower_left_m,
                total_rr_m3,
                activation_dataset,
            )

            self.activation_file = activation_file_path

        ## common workflow for both types of simulations

        # apply rototranslation to the reaction rate mesh
        if (
            not np.array_equal(arg_dict["activation_rotation_degs"], [0.0, 0.0, 0.0])
        ) or (
            not np.array_equal(arg_dict["activation_translation_m"], [0.0, 0.0, 0.0])
        ):
            if not self.activation_file:
                raise ValueError(
                    "cannot apply rototranslation to activation dataset if \
                    activation_file is not specified"
                )

            print("applying rototranslation to activation dataset ...")

            self.activation_rotation_center_mode = arg_dict[
                "activation_rotation_center_mode"
            ]
            self.activation_rotation_euler_order = arg_dict[
                "activation_rotation_euler_order"
            ]
            self.activation_rotation_degs = arg_dict["activation_rotation_degs"]
            self.activation_translation_m = arg_dict["activation_translation_m"]

            input_vtk = Path(self.activation_file)
            output_vtk = input_vtk.with_name(f"{input_vtk.stem}_rototrans.vtk")

            apply_roto_translation_to_vtk_grid(
                input_vtk,
                output_vtk,
                self.activation_rotation_degs,
                self.activation_translation_m,
                order=self.activation_rotation_euler_order,
                center=self.activation_rotation_center_mode,
            )

            self.activation_file = output_vtk

        if "activation_constant" in arg_dict:
            self.activation_constant = float(arg_dict["activation_constant"])

        if "activation_normalization" in arg_dict:
            self.activation_normalization = float(arg_dict["activation_normalization"])

        if "fv_scheme" in arg_dict:
            if arg_dict["fv_scheme"].lower() in ["stable", "accurate"]:
                self.fv_scheme = arg_dict["fv_scheme"]
            else:
                raise ValueError("fv_scheme must be 'stable' or 'accurate'")

        # set cfd type
        cfd_type = arg_dict.get("cfd_type")
        if cfd_type is None:
            raise ValueError("ERROR: cfd type not specified in input file")
        self.cfd_type = cfd_type.lower()

        # set fluent fluid region name
        if self.cfd_type == "fluent-h5-multi":
            region = arg_dict.get("fluent_fluid_region_name")
            if region is None:
                raise ValueError(
                    "ERROR: name of the fluid region to extract not specified! "
                    "Use parameter FLUENT_FLUID_REGION_NAME"
                )
            self.fluent_fluid_region_name = region

        # set time treatment
        time_treatment = arg_dict.get("time_treatment")
        if time_treatment is None:
            raise ValueError("ERROR: type of time treatment not specified")
        if time_treatment.lower() not in {"steadystate", "transient"}:
            raise ValueError(
                "ERROR: time treatment must be 'steadyState' or 'transient'"
            )
        self.time_treatment = time_treatment.lower()

        # set transport parameters
        self.molecular_diffusion = float(arg_dict["molecular_diffusion"])
        self.schmidt_number = float(arg_dict["schmidt_number"])
        self.inlet_conc = float(arg_dict["inlet_conc"])

        # set cfd and fluned paths
        self.cfd_path = Path(arg_dict["cfd_path"])
        if not self.cfd_path.is_dir():
            raise FileNotFoundError(f"Folder not found: {self.cfd_path}")

        self.fluned_path = Path(self.cfd_path) / self.case

        return

    def generate_fluned_files(self):
        """
        this function initialize the cfd class, after the initial openfoam or
        fluent simulation has been processed
        """

        self.cfd_simulation.post_process_openfoam_simulation()

        self.fluned_simulation = oFoamFluned(self.fluned_path)

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
            self.activation_constant,
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
        it initialize a oFoamBase object and parse its data
        """

        self.fluned_simulation = oFoamFluned(self.fluned_path)

        self.fluned_simulation.post_process_openfoam_simulation()
        self.fluned_simulation.post_process_fluned_simulation()

        return

    def initialize_cfd_class(self):
        """
        this function initialize the cfd class, after the initial openfoam or
        fluent simulation has been processed
        """

        self.cfd_simulation = oFoamBase(self.cfd_path)

        return

    def create_case_folders_fluent(self):
        """
        This function creates additional folders required for the FLUENT simulation converted
        into the OpenFOAM format
        """

        case_folder = self.cfd_path

        # T=1 folder
        one_folder = Path(case_folder) / "1"
        one_folder.mkdir(exist_ok=True)

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
        target_file = Path(self.fluned_path) / "0" / "phi"
        origin_file = Path(self.cfd_path) / str(last_time) / "phi"

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
        target_file = Path(self.fluned_path) / "0" / "U"
        origin_file = Path(self.cfd_path) / str(last_time) / "U"
        shutil.copyfile(origin_file, target_file)

        return None

    def copy_last_nut(self):
        """this function look for the last nut file in the cfd folder,
        if it does not exists it means the simulation is laminar"""

        last_time = self.cfd_simulation.last_time

        targetFile = Path(self.fluned_path) / "0" / "nut"
        originFile = Path(self.cfd_path) / str(last_time) / "nut"
        checkFile = originFile.is_file()
        if checkFile:
            shutil.copyfile(originFile, targetFile)

        return

    def copy_poly_mesh(self):
        """this function copy the poly Mesh files in the FLUNED folder"""

        sourceFolder = Path(self.cfd_path) / "constant" / "polyMesh"
        targetFolder = Path(self.fluned_path) / "constant" / "polyMesh"

        # shutil do not allow to copy into an existing folder
        dirCheck = targetFolder.is_dir()
        if dirCheck:
            shutil.rmtree(targetFolder)

        shutil.copytree(sourceFolder, targetFolder)

        return

    def launch_check_mode(self):
        """
        this function launches some functional objects and read its results for some additional checks
        """

        launch_grad_func_object(self.fluned_path)
        self.fluned_simulation.read_velocities()

        return

    def write_results(self, arguments):
        """
        this function generate the final results
        """

        if self.fluned_simulation.time_treatment == "steadystate":
            self.write_summary_steady(arguments)
            self.write_summary_xml()
        else:
            self.write_summary_transient()

        return None

    def write_summary_transient(self):
        """
        This function writes the summary file in the RESULTS/ folder
        """

        print(
            "WARNING: summary writing for transient simulations is not implemented"
            "due to ongoing refactoring"
            "contact the developer for more info"
        )

        return

    def write_summary_steady(self, arguments):
        """
        This function writes the summary file in the RESULTS/ folder
        """

        summary_file = Path(self.fluned_simulation.results_folder) / "SUMMARY.csv"

        results_sim = self.fluned_simulation

        inlet_atoms = abs(results_sim.total_inlet_t_atoms)
        inlet_activity = abs(results_sim.inlet_td_bq_m3[-1])

        outlet_activity = results_sim.outlet_t_bq_m3[-1]
        outlet_atoms = results_sim.total_outlet_t_atoms

        tot_activity = results_sim.total_isotope_activity
        avg_activity = results_sim.total_average_isotope_activity

        averageVolume = sum(results_sim.volumes) / len(results_sim.volumes)

        faces_post = sorted(results_sim.patches.values(), key=lambda x: x.face_id)

        if arguments.check:
            decay_const = self.decay_constant
            tot_scalar = sum(results_sim.t_scalar)

            # parameter 2 based on transit times
            transit_times = [
                (vol ** (1 / 3)) / mod(v)
                for vol, v in zip(results_sim.volumes, results_sim.velocities)
            ]
            average_transit_time = sum(transit_times) / len(transit_times)

            mesh_quality_parameter_2 = [
                t * c * decay_const / (np.log(2) * tot_scalar)
                for t, c in zip(transit_times, results_sim.t_scalar)
            ]
            average_mesh_quality_parameter_2 = sum(mesh_quality_parameter_2) / len(
                mesh_quality_parameter_2
            )

        negCheck = any(n < 0 for n in results_sim.t_scalar)

        with open(summary_file, "w") as fw:
            fw.write("FLUNED SIMULATION SUMMARY\n")
            fw.write("CASE,{},\n".format(results_sim.case))
            fw.write("STEADY STATE SIMULATION\n")
            fw.write("N ELEMENTS,{},\n".format(results_sim.n_internal_cells))
            fw.write("ISOTOPE,{},\n".format(results_sim.isotope.upper()))
            fw.write("DECAY CONSTANT,{:e},\n".format(results_sim.decay_constant))
            fw.write("CASE VOLUME [m3],{:e},\n".format(results_sim.volume_m3))
            fw.write("MOL DIFFUSION,{:e},\n".format(results_sim.molecular_diffusion))
            fw.write("TURB SCHMIDT N,{:f},\n".format(results_sim.schmidt_number))
            fw.write("\n")
            fw.write("\n")
            fw.write("QUALITY\n")
            fw.write("AVG VOL [m3],{:e},\n".format(averageVolume))
            if arguments.check:
                fw.write("AVG CELL TRNST T [s],{:f},\n".format(average_transit_time))  # type: ignore
                fw.write("AVG CELL Q2,{:e},\n".format(average_mesh_quality_parameter_2))  # type: ignore
            if negCheck:
                fw.write("WARNING some elements are negative")
            if arguments.cdgs:
                fw.write("\n")
                fw.write("\n")
                fw.write("source sampling\n")
                fw.write(
                    "sampling resolution [m],{:f},\n".format(
                        self.source_sampling_resolution_cm / 100
                    )
                )
                fw.write(
                    "sampled voxels [#],{:d},\n".format(
                        results_sim.x_ints * results_sim.y_ints * results_sim.z_ints
                    )
                )
                fw.write(
                    "vtk emission rate [#/s],{:e},\n".format(
                        results_sim.total_isotope_emission_rate
                    )
                )
                sampstring = "sampled emission rate (unscaled) [#/s],{:e},\n"
                fw.write(sampstring.format(results_sim.raw_sampled_total_emission_rate))
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
                string = "INLET-NORMALIZED AVG ACTIVITY,{:.5e},\n"
                fw.write(string.format(norm_avg_activity))

            if inlet_atoms != 0:
                reduction_rate = abs(outlet_atoms / inlet_atoms)
                fw.write("OUT/IN RATIO,{:.5e},\n".format(reduction_rate))

            fw.write("\n")
            fw.write("\n")
            fw.write("FACES\n")

            for face in faces_post:
                if face.face_type != "wall":
                    fw.write(face.to_csv_patch_element())

        return

    def write_summary_xml(self):
        """
        Write a compact XML summary file in the RESULTS/ folder
        with selected simulation metadata, quality, and activation data.
        """

        xml_file = Path(self.fluned_simulation.results_folder) / "fluned_summary.xml"
        results_sim = self.fluned_simulation

        inlet_atoms = abs(results_sim.total_inlet_t_atoms)
        inlet_activity = abs(results_sim.inlet_td_bq_m3[-1])

        outlet_atoms = results_sim.total_outlet_t_atoms
        outlet_activity = results_sim.outlet_t_bq_m3[-1]

        tot_activity = results_sim.total_isotope_activity
        avg_activity = results_sim.total_average_isotope_activity

        average_volume = sum(results_sim.volumes) / len(results_sim.volumes)

        faces_post = sorted(results_sim.patches.values(), key=lambda x: x.face_id)

        root = ET.Element("simulation_summary")
        root.set("simulation_type", "steady_state")

        ET.SubElement(root, "case_name").text = str(results_sim.case)
        ET.SubElement(root, "number_of_elements").text = str(
            results_sim.n_internal_cells
        )
        ET.SubElement(root, "isotope").text = str(results_sim.isotope.upper())

        volume = ET.SubElement(root, "mesh_volume", unit="m^3")
        volume.text = f"{results_sim.volume_m3:.6e}"

        decay = ET.SubElement(root, "decay_constant", unit="s^-1")
        decay.text = f"{results_sim.decay_constant:.6e}"

        diffusion = ET.SubElement(
            root, "molecular_diffusion_coefficient", unit="m^2 s^-1"
        )
        diffusion.text = f"{results_sim.molecular_diffusion:.6e}"

        schmidt = ET.SubElement(root, "turbulent_schmidt_number")
        schmidt.text = f"{results_sim.schmidt_number:.6f}"

        # QUALITY section
        quality = ET.SubElement(root, "quality")
        avg_vol = ET.SubElement(quality, "average_cell_volume", unit="m^3")
        avg_vol.text = f"{average_volume:.6e}"

        # ACTIVATION section
        activation = ET.SubElement(root, "activation")

        ET.SubElement(
            activation, "inlet_atoms", unit="# s^-1"
        ).text = f"{inlet_atoms:.5e}"
        ET.SubElement(
            activation, "outlet_atoms", unit="# s^-1"
        ).text = f"{outlet_atoms:.5e}"
        ET.SubElement(
            activation, "total_inlet_activity", unit="Bq m^-3"
        ).text = f"{inlet_activity:.5e}"
        ET.SubElement(
            activation, "total_outlet_activity", unit="Bq m^-3"
        ).text = f"{outlet_activity:.5e}"
        ET.SubElement(
            activation, "total_case_activity", unit="Bq"
        ).text = f"{tot_activity:.5e}"
        ET.SubElement(
            activation, "total_average_activity", unit="Bq m^-3"
        ).text = f"{avg_activity:.5e}"

        if inlet_activity != 0:
            norm_avg_activity = abs(avg_activity / inlet_activity)
            ET.SubElement(
                activation, "inlet_normalized_average_activity"
            ).text = f"{norm_avg_activity:.5e}"

        if inlet_atoms != 0:
            out_in_ratio = abs(outlet_atoms / inlet_atoms)
            ET.SubElement(
                activation, "outlet_to_inlet_ratio"
            ).text = f"{out_in_ratio:.5e}"

        # PATCH section
        patches = ET.SubElement(root, "patches")

        for patch in faces_post:
            if patch.face_type != "wall":
                patches.append(patch.to_xml_patch_element())

        tree = ET.ElementTree(root)
        ET.indent(tree, space="  ", level=0)
        tree.write(xml_file, encoding="utf-8", xml_declaration=True)

        return
