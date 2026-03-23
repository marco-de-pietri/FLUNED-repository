import math
import re
import shutil
import warnings
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pyvista as pv

from fluned.isotopes.isotopes import bins_from_lines, load_isotopes

from .fluned_tool_launchers import launch_centroid_func_object, launch_foam_to_vtk
from .fluned_vtk_utils import (
    generate_external_stl,
    generate_triangularized_h5m_um_mesh,
    generate_triangularized_scalar_mesh,
    get_vtk_celldata_array,
    get_vtk_dimensions,
    get_vtk_volumes,
    sample_vtk,
    write_cartesian_vtk,
)
from .ofoam_base import oFoamBase
from .patch_class import get_post_process_list


def formatValues(vector):
    maxLen = 70
    returnString = ""
    newLine = ""

    for item in vector:
        newNumber = "{:.7e}".format(item)
        if returnString == "" and newLine == "":
            newLine = newNumber
            continue

        if len(newLine + " " + newNumber) > maxLen:
            returnString += newLine + "\n"
            newLine = newNumber

        else:
            newLine = newLine + " " + newNumber

    returnString += newLine + "\n"

    return returnString


def bounding_box_corners(points, padding=1.0):
    """
    Given a list of coordinates (each a list/tuple of equal length),
    return the two opposite corners of the minimal axis-aligned
    bounding box that encompasses all points, with optional padding.

    Args:
        points: List of coordinates (each a list/tuple of equal length)
        padding: Padding to add to all sides of the bounding box.
                Can be a single number (applied to all dimensions) or
                a list/tuple of numbers (one per dimension).

    Returns:
        (min_corner, max_corner), where each is a tuple of length D.
    """
    if not points:
        raise ValueError("points must be a non-empty list")

    dim = len(points[0])
    for p in points:
        if len(p) != dim:
            raise ValueError("all points must have the same dimension")

    # Handle padding parameter
    if isinstance(padding, (int, float)):
        # Single padding value for all dimensions
        padding_values = [padding] * dim
    else:
        # Padding per dimension
        padding_values = list(padding)
        if len(padding_values) != dim:
            raise ValueError(
                "padding must be a single number or have same length as point dimensions"
            )

    mins = [float("inf")] * dim
    maxs = [float("-inf")] * dim

    for p in points:
        for i, v in enumerate(p):
            if v < mins[i]:
                mins[i] = v
            if v > maxs[i]:
                maxs[i] = v

    # Calculate dimensions and apply padding as scaling factor
    padded_mins = []
    padded_maxs = []

    for i in range(dim):
        dimension_size = maxs[i] - mins[i]
        padding_amount = dimension_size * padding_values[i]

        # Handle edge case where dimension_size is 0 (all points identical in this dimension)
        if dimension_size == 0:
            padding_amount = abs(
                padding_values[i]
            )  # Use absolute padding for zero-width dimensions

        padded_mins.append(mins[i] - padding_amount)
        padded_maxs.append(maxs[i] + padding_amount)

    return tuple(padded_mins), tuple(padded_maxs)


def mesh_ints_from_bounds(
    lower: Iterable[float],
    upper: Iterable[float],
    min_voxel_size: float | Sequence[float],
    max_voxels_n: int | None = None,
) -> tuple[int, ...]:
    """
    Compute the number of cells (intervals) per axis for a Cartesian mesh so that
    the actual voxel size along each axis is <= the provided minimum, and (optionally)
    the total number of voxels is capped by scaling the intervals down.

    Parameters
    ----------
    lower : iterable of numbers
        Lower-left (min) coordinates per dimension, e.g. (x0, y0, z0).
    upper : iterable of numbers
        Upper-right (max) coordinates per dimension, e.g. (x1, y1, z1).
    min_voxel_size : number or sequence of numbers
        Minimum voxel size. If a single number, it's applied to all dimensions.
        If a sequence, it must have the same length as `lower`/`upper`.
    max_voxels_n : int or None, optional
        Maximum allowed total number of voxels (product of intervals across axes).
        If the unconstrained mesh would exceed this cap, intervals are scaled down
        approximately uniformly across axes (preserving aspect ratio as much as possible)
        until the product is <= max_voxels_n.

    Returns
    -------
    tuple of int
        Number of intervals (cells) along each dimension.

    Notes
    -----
    - Base counts use ceil((upper - lower) / min_size) per dimension, with a minimum of 1.
    - If (upper <= lower) along any axis, a ValueError is raised.
    - If max_voxels_n is set and the base total exceeds it, voxel sizes may end up
      larger than min_voxel_size due to the global cap.
    """
    lo = tuple(float(v) for v in lower)
    hi = tuple(float(v) for v in upper)
    if len(lo) != len(hi):
        raise ValueError("`lower` and `upper` must have the same length.")

    dim = len(lo)

    if isinstance(min_voxel_size, (int, float)):
        mins = (float(min_voxel_size),) * dim
    else:
        mins = tuple(float(v) for v in min_voxel_size)
        if len(mins) != dim:
            raise ValueError(
                "min_voxel_size must be a scalar or a sequence with the same length as lower/`upper"
            )

    base_counts: list[int] = []
    for a, b, m in zip(lo, hi, mins, strict=True):
        if m <= 0:
            raise ValueError("All minimum voxel sizes must be > 0.")
        if b <= a:
            raise ValueError(
                "Each upper bound must be greater than the corresponding lower bound."
            )
        n = max(1, math.ceil((b - a) / m))
        base_counts.append(int(n))

    # No cap requested
    if max_voxels_n is None:
        return tuple(base_counts)

    if not isinstance(max_voxels_n, int) or max_voxels_n < 1:
        raise ValueError("`max_voxels_n` must be an int >= 1 or None.")

    def prod(xs: Sequence[int]) -> int:
        p = 1
        for v in xs:
            p *= int(v)
        return p

    base_total = prod(base_counts)
    if base_total <= max_voxels_n:
        return tuple(base_counts)

    # Uniform scale factor in "interval-count space" to hit the volume cap
    scale = (max_voxels_n / base_total) ** (1.0 / dim)

    # Start with floored scaled counts to ensure we don't overshoot too easily
    counts = [max(1, int(math.floor(c * scale))) for c in base_counts]

    # In rare cases (e.g., scale ~1 and flooring didn't reduce enough), ensure product <= cap.
    # Decrease counts one-by-one, prioritizing the axes with the largest current counts.
    while prod(counts) > max_voxels_n:
        # pick axis with largest count that can still be reduced
        k = max(
            (i for i, v in enumerate(counts) if v > 1),
            key=lambda i: counts[i],
            default=None,
        )
        if k is None:
            break  # cannot reduce further
        counts[k] -= 1

    return tuple(counts)


class oFoamFluned(oFoamBase):
    def __init__(self, path: str | Path):
        super().__init__(path)
        self.reduction_rate = []
        self.normalized_average_td = []
        self.inlet_td_conc_atoms_m3 = []
        self.outlet_rr_conc_atoms_m3 = []
        self.average_rr_conc_atoms_m3 = []
        self.post_process_td_average = []
        self.average_ta = []
        self.volume_m3 = 0
        self.vtk_file_path = ""
        self.scaled_vtk_file_path = ""
        self.vtk_dimensions = []

    def post_process_fluned_simulation(self):
        """
        this function is called when we process a finished FLUNED simulation
        """
        for face in self.patches.values():
            face.post_process_face()

        self.post_process_td_average, _, _ = get_post_process_list(
            self.post_process_path, "volTdSum", "", "volAverage(Td)"
        )
        self.average_ta, _, _ = get_post_process_list(
            self.post_process_path, "volTaSum", "", "volAverage(Ta)"
        )

        self.total_inlet_t_atoms = self.get_total_inlet_t_atoms()
        self.total_outlet_t_atoms = self.get_total_outlet_t_atoms()

        self.inlet_td_conc_atoms_m3 = self.get_inlet_td_conc_atoms_m3()
        self.outlet_t_conc_atoms_m3 = self.get_outlet_t_conc_atoms_m3()
        # print ("self.inlet_td_conc_atoms_m3")
        # print (f"{self.inlet_td_conc_atoms_m3[-1]:.2e}")

        self.reduction_rate = self.get_reduction_rate()
        # print ("self.reduction_rate")
        # print (self.reduction_rate)

        self.normalized_average_td = self.get_normalized_average_td()
        # print ("self.normalized_average_td")
        # print (self.normalized_average_td)

        self.outlet_rr_conc_atoms_m3 = self.get_outlet_rr_conc_atoms_m3()
        # print ("self.outlet_rr_conc_atoms_m3")
        # print (f"{self.outlet_rr_conc_atoms_m3[-1]:.2e}")

        # print ("self.average_ta")
        # print (f"{self.average_ta[-1]:.2e}")

        self.parse_constants_file()
        self.assign_isotope_data()
        self.get_time_treatment()
        self.read_volumes()
        self.read_t()

        self.results_folder = self.path / "RESULTS"
        self.results_folder.mkdir(exist_ok=True)

        launch_foam_to_vtk(str(self.path))
        self.vtk_file_folder = self.path / "VTK"
        self.get_vtk_file()
        self.vtk_dimensions, self.volume_m3 = get_vtk_dimensions(self.vtk_file_path)

        generate_external_stl(self.vtk_file_path, self.results_folder)

        self.isotope_amounts = [
            vol * t for vol, t in zip(self.volumes, self.t_scalar)
        ]  # atoms/cell
        self.total_isotope_amount = sum(self.isotope_amounts)  # atoms

        self.total_average_isotope_concentration = (
            self.total_isotope_amount / self.volume_m3
        )
        self.total_isotope_activity = (
            self.decay_constant * self.total_isotope_amount
        )  # decays/s
        self.total_isotope_emission_rate = (
            self.tot_p_emission * self.total_isotope_activity
        )  # decay particle/s

        return None

    def get_vtk_file(self, custom_vtk: str | None = None):
        """
        this function looks for the vtk file in the VTK folder generated by the
        vtk openfoam utility and store it for generation of results

        custom values can be provided for special workflow where the vtk
        simulation file is rototranslated before generating the radiation source
        """

        vtkFiles = []

        self.vtk_path = ""

        if custom_vtk is not None:
            filename = Path(custom_vtk).stem
        else:
            filename = ""

        pat_string = r"{}\.vtk\Z".format(filename)

        vtkFilePat = re.compile(pat_string, re.IGNORECASE)

        for vtk_path in self.vtk_file_folder.iterdir():
            filename = vtk_path.name
            matchVTKfiles = vtkFilePat.findall(filename)
            if len(matchVTKfiles) == 1:
                vtkFiles.append(filename)

        if len(vtkFiles) == 1:
            self.vtk_file_path = self.path / "VTK" / vtkFiles[0]
        else:
            print("ERROR zero or more than one vtk files")
            print("found in the VTK folder")
            print("found vtk files: ", vtkFiles)
            raise ValueError()

        return

    def get_reduction_rate(self):
        """
        this function calculates the reduction rate of the isotope in the cfd
        it sums all the outlet Td sum values - therefore multiple inlets or
        outlets are not supported yet
        """

        atom_in = [0.0]
        atom_out = [0.0]

        for patch in self.patches.values():
            if patch.face_type == "inlet":
                atom_in = patch.post_process_td_flow
            if patch.face_type == "outlet":
                atom_out = patch.post_process_td_flow

        red_ratio = [abs(x / y) if y != 0 else 0 for x, y in zip(atom_out, atom_in)]

        return red_ratio

    def get_outlet_rr_conc_atoms_m3(self):
        """
        this function returns the concentration of the outlet patch for the
        Ta field, meaning for the part generated only by the irradiation and
        not from the inlet flow
        at the moment it supports only one outlet patch
        """

        conc_out = []

        for patch in self.patches.values():
            if patch.face_type == "outlet":
                conc_out = patch.ta_conc_atoms_m3
                break

        return conc_out

    def get_inlet_td_conc_atoms_m3(self):
        """
        this function returns the concentration of the inlet patch for the
        Td field, meaning for the part due to the inlet flow

        at the moment it supports only one outlet patch
        """

        conc_out = []

        for patch in self.patches.values():
            if patch.face_type == "inlet":
                conc_out = patch.td_conc_atoms_m3
                break

        return conc_out

    def get_outlet_t_conc_atoms_m3(self):
        """
        this function returns the concentration of the outlet patch for the
        T field

        at the moment it supports only one outlet patch
        """

        conc_out = []

        for patch in self.patches.values():
            if patch.face_type == "outlet":
                conc_out = patch.t_conc_atoms_m3
                break

        return conc_out

    def get_total_inlet_t_atoms(self):
        """
        this function returns the total amount of atoms
        enetering the mesh through the inlets
        """

        atoms_in = 0

        for patch in self.patches.values():
            if patch.face_type == "inlet":
                atoms_in += patch.post_process_t_flow[-1]

        return atoms_in

    def get_total_outlet_t_atoms(self):
        """
        this function returns the total amount of atoms
        exiting the mesh through the inlets
        """

        atoms_out = 0

        for patch in self.patches.values():
            if patch.face_type == "outlet":
                atoms_out += patch.post_process_t_flow[-1]

        return atoms_out

    def get_normalized_average_td(self):
        """
        this function returns the average concentration of the Td field due
        only to the decay of the inlet flow
        """

        norm_td = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_td_average, self.inlet_td_conc_atoms_m3)
        ]

        return norm_td

    def run_cartesian_sampling(
        self,
        source_sampling_resolution_cm: float | None = None,
        sampling_dataset: str | None = None,
    ):
        """
        this function sample the fluned results into a cartesian mesh
        """

        self.calculate_cartesian_sampling_coordinates(source_sampling_resolution_cm)
        self.sample_source_to_cartesian_mesh(sampling_dataset)

        self.write_sampled_cartesian_source_vtk()

        return

    def scale_mesh_results(self, out_path, inlet_activity, decay_const):
        """
        this function takes the vtk um-mesh as computed by fluned, it scales
        it to a new inlet activity, and writes the scaled results to a new vtk file

        IMPORTANT: it converts concentrations to activities
        """

        print("writing vtk file for circuit calculation ... ")

        self.scaled_vtk_file_path = out_path

        mesh = pv.read(self.vtk_file_path)

        decay_array = mesh.cell_data["Td"]
        rr_array = mesh.cell_data["Ta"]

        del mesh.cell_data["Td"]
        del mesh.cell_data["Ta"]
        del mesh.cell_data["T"]
        # del mesh.cell_data['CellID']
        del mesh.point_data["T"]
        del mesh.point_data["Ta"]
        del mesh.point_data["Td"]

        decay_array = inlet_activity * decay_array / self.inlet_td_conc_atoms_m3[-1]
        rr_array = decay_const * rr_array

        mesh.cell_data["average_vol_activity_bq_m3_decay"] = decay_array
        mesh.cell_data["average_vol_activity_bq_m3_rr"] = rr_array
        mesh.cell_data["average_vol_activity_bq_m3"] = rr_array + decay_array

        mesh.save(out_path)

        return

    def calculate_cartesian_sampling_coordinates(self, sampling_res_cm: float | None):
        """
        this function returns a list of dictionaries with the info relative to
        the sampling coordinates
        """

        if sampling_res_cm is None or sampling_res_cm <= 0:
            raise ValueError("sampling_res_cm argument must be a positive number")

        vtk_boundaries = self.vtk_dimensions

        # convert to cm
        x_bounds = [
            math.floor(vtk_boundaries[0] * 100),
            math.ceil(vtk_boundaries[1] * 100),
        ]
        y_bounds = [
            math.floor(vtk_boundaries[2] * 100),
            math.ceil(vtk_boundaries[3] * 100),
        ]
        z_bounds = [
            math.floor(vtk_boundaries[4] * 100),
            math.ceil(vtk_boundaries[5] * 100),
        ]

        x_ints = math.ceil((x_bounds[1] - x_bounds[0]) / (sampling_res_cm))
        y_ints = math.ceil((y_bounds[1] - y_bounds[0]) / (sampling_res_cm))
        z_ints = math.ceil((z_bounds[1] - z_bounds[0]) / (sampling_res_cm))

        x_nodes = [(x_bounds[0] + sampling_res_cm * i) for i in range(x_ints + 1)]
        y_nodes = [(y_bounds[0] + sampling_res_cm * i) for i in range(y_ints + 1)]
        z_nodes = [(z_bounds[0] + sampling_res_cm * i) for i in range(z_ints + 1)]

        sample_coordinates_m = []

        voxel_list = []
        voxel_id = 1

        for i in range(x_ints):
            # xVoxelNodes = [x_nodes[i], x_nodes[i+1]]
            x_voxel_center = (x_nodes[i + 1] + x_nodes[i]) / 2
            for j in range(y_ints):
                # yVoxelNodes = [y_nodes[j], y_nodes[j+1]]
                y_voxel_center = (y_nodes[j + 1] + y_nodes[j]) / 2
                for k in range(z_ints):
                    # zVoxelNodes = [z_nodes[k], z_nodes[k+1]]
                    z_voxel_center = (z_nodes[k + 1] + z_nodes[k]) / 2
                    new_dict = {}
                    new_dict["id"] = voxel_id
                    coords = [x_voxel_center, y_voxel_center, z_voxel_center]
                    new_dict["cent_coords_cm"] = coords
                    coords_m = [coord / 100 for coord in coords]
                    new_dict["cent_coords_m"] = coords_m
                    sample_coordinates_m.append(coords_m)
                    new_dict["volume_cm3"] = sampling_res_cm**3
                    voxel_list.append(new_dict)
                    voxel_id += 1

        self.x_ints = x_ints
        self.y_ints = y_ints
        self.z_ints = z_ints

        self.x_nodes = x_nodes
        self.y_nodes = y_nodes
        self.z_nodes = z_nodes
        self.cartesian_voxel_list = voxel_list
        self.cartesian_sample_coordinates_m = sample_coordinates_m
        self.cartesian_sampling_precision_cm = sampling_res_cm
        self.cartesian_voxel_volume_cm3 = sampling_res_cm**3

        return

    def sample_source_to_cartesian_mesh(self, dataset: str | None):
        """
        this function takes the cartesian sampling coordinates computed before and
        use them to:
        1. sample the concentrations from the fluned simulation um mesh
        2. calculate the total sampled values
        3. the sampled emission rates are adjusted so the total emission remains the same

        NB the generated sampled values are in cm scale as this process is done to
        generate radiation sources either for openMC or MCNP
        """

        if dataset is None or dataset == "":
            raise ValueError("dataset argument cannot be None or empty string")

        sampled_cartesian_concentrations_cm3 = sample_vtk(
            self.vtk_file_path, dataset, self.cartesian_sample_coordinates_m
        )

        voxel_volume = self.cartesian_sampling_precision_cm**3

        for voxel, concentration in zip(
            self.cartesian_voxel_list, sampled_cartesian_concentrations_cm3
        ):
            voxel["emission_rate_per_voxel"] = (
                concentration
                * voxel_volume
                * self.decay_constant
                * self.tot_p_emission
                * 1e-06
            )  # atoms per m3 to cm3

        self.raw_sampled_total_emission_rate = sum(
            [voxel["emission_rate_per_voxel"] for voxel in self.cartesian_voxel_list]
        )

        ratio_vtk_sampling = (
            self.total_isotope_emission_rate / self.raw_sampled_total_emission_rate
        )

        for voxel in self.cartesian_voxel_list:
            voxel["normalized_emission_rate_per_voxel"] = (
                voxel["emission_rate_per_voxel"] * ratio_vtk_sampling
            )

        normalization_check = sum(
            [
                voxel["normalized_emission_rate_per_voxel"]
                for voxel in self.cartesian_voxel_list
            ]
        )

        if not math.isclose(
            self.total_isotope_emission_rate,
            normalization_check,
            rel_tol=1e-03,
            abs_tol=0.0,
        ):
            print("Error in emission rate sampling normalization")
            print("total isotope emission rate: ", self.total_isotope_emission_rate)
            print("sum of normalized emission rates: ", normalization_check)
            raise ValueError()

        return

    def generate_zero_t(self, inlet_conc):
        """
        this function generate the T file at t=0
        """

        zero_folder = self.path / "0"
        zero_t_path = zero_folder / "T"

        t_header_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      T;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}

"""

        with open(zero_t_path, "w", encoding="utf-8") as fw:
            fw.write(t_header_text)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                if face.face_type == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "inlet":
                    fw.write("        type            fixedValue;\n")
                    vals_string = "        value         uniform {};\n"
                    fw.write(vals_string.format(inlet_conc))
                else:
                    raise ValueError(
                        "ERROR: patch type not recognized - ", face.face_type
                    )

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_zero_ta(self):
        """
        this function generate the Ta file at t=0
        """

        zero_folder = self.path / "0"
        zero_ta_path = zero_folder / "Ta"

        ta_header_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Ta;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}
"""

        with open(zero_ta_path, "w", encoding="utf-8") as fw:
            fw.write(ta_header_text)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                if face.face_type == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "inlet":
                    fw.write("        type            fixedValue;\n")
                    vals_string = "        value          uniform 0;\n"
                    fw.write(vals_string)
                else:
                    raise ValueError(
                        "ERROR: patch type not recognized - ", face.face_type
                    )

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_zero_tr(self, time_treatment):
        """
        this function generate the Tr file at t=0
        """

        zero_folder = self.path / "0"
        zero_tr_path = zero_folder / "Tr"

        tr_header_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Tr;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}
"""

        with open(zero_tr_path, "w", encoding="utf-8") as fw:
            fw.write(tr_header_text)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                if face.face_type == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "inlet":
                    fw.write("        type            fixedValue;\n")
                    if time_treatment == "steadystate":
                        vals_string = "        value          uniform 0;\n"
                    elif time_treatment == "transient":
                        vals_string = "        value          uniform 1;\n"
                    else:
                        raise ValueError("time_treatment argument not recognized")
                    fw.write(vals_string)
                else:
                    raise ValueError(
                        "ERROR: patch type not recognized - ", face.face_type
                    )

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_zero_td(self, inlet_conc):
        """
        this function generate the Td file at t=0
        """

        zero_folder = self.path / "0"
        zero_td_path = zero_folder / "Td"

        td_header_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Td;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}

"""

        with open(zero_td_path, "w", encoding="utf-8") as fw:
            fw.write(td_header_text)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                if face.face_type == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "inlet":
                    fw.write("        type            fixedValue;\n")
                    vals_string = "        value         uniform {};\n"
                    fw.write(vals_string.format(inlet_conc))
                else:
                    raise ValueError(
                        "ERROR: patch type not recognized - ", face.face_type
                    )

                fw.write("    }\n")

            fw.write(end_text)

        return

    def assign_activation_rates(
        self,
        activation_file,
        activation_constant,
        activation_dataset,
        activation_dataset_error,
        activation_normalization,
    ):
        """this file create computes the reaction rates that are later written in the Source
        file.
        NB the source rates are intended as volumetric [#.m-3.s-1]

        if the simulation was done with openmc or mcnp with no scale applied
        it must be taken into account that the reaction rates are computed over
        a mesh with volumetric values computed in cm [#.cm-3.s-1]

        Fluned does not apply any correction factor for this, either apply a scaling factor in
        the fluned input or modify the input data accordingly

        There are three modes:
        1) no activation, then the file contains only zeros.
        2) constant activation (with input value)
        3) with a source file
        """

        if activation_file == "" or activation_file is None:
            if activation_constant is None or math.isclose(
                activation_constant, 0.0, abs_tol=1e-12
            ):
                # case with no rad source
                activ_sources = [0 for _ in range(self.n_internal_cells)]
            else:
                # case with constant value source
                self.read_volumes()
                activ_sources = [activation_constant * vol for vol in self.volumes]

        else:
            activation_file = Path(activation_file)
            # case with source file
            launch_centroid_func_object(self.path)
            self.read_volumes()
            self.read_centroids()

            # 1.sample the activation file
            print("Sampling Reaction Rate file ... ")
            sampledRates = sample_vtk(
                activation_file, activation_dataset, self.centroids
            )
            # remove possible small negative values
            sampledRates = [val if val > 0 else 0 for val in sampledRates]

            # 1.1 if present sample the vtk file to get the error array
            if activation_dataset_error != "":
                print("Sampling Reaction Rate MCNP errors  ... ")
                self.sampled_rr_error = sample_vtk(
                    activation_file, activation_dataset_error, self.centroids
                )

            # 2.use the activation const as a factor
            if float(activation_constant) == 0.0:
                factor = 1
            else:
                factor = activation_constant
            activ_sources = [factor * rate for rate in sampledRates]

            # print("debug: total sampled reaction rate #/s pre scaling")
            # print(sum([rate * vol for rate, vol in zip(sampledRates, self.volumes)]))

            # print("debug: total sampled reaction rate #/s after scaling")
            # print(sum([rate * vol for rate, vol in zip(activ_sources, self.volumes)]))

            # 3.apply a normalization factor if provided
            if activation_normalization != 0:
                vec = [rate * vol for rate, vol in zip(activ_sources, self.volumes)]

                total_sampled = sum(vec)
                print("total sampled atoms/s")
                print(total_sampled)
                if total_sampled == 0:
                    raise ValueError(
                        "ERROR: total sampled reaction rate is zero, cannot apply normalization factor"
                    )
                normalization_factor = activation_normalization / total_sampled

                activ_sources = [rate * normalization_factor for rate in activ_sources]

                nVec = [rate * vol for rate, vol in zip(activ_sources, self.volumes)]

                print("new total sampled atoms/s")
                print(sum(nVec))

        self.reaction_rates = activ_sources

        return None

    def generate_source_file(self, activation_dataset_error=""):
        """
        generate the source file using the dataset reaction rate data
        """

        zeroFolder = self.path / "0"
        zeroSourcePath = zeroFolder / "Source"

        sHeaderText = """
        FoamFile
        {
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0";
            object      Source;
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        dimensions      [0 0 -1 0 0 0 0];

        internalField   nonuniform List<scalar>
        """
        boundaryText = """
        boundaryField
        {

        """
        closerText = """

        }


        // *********************************************************************** //
        """

        with open(zeroSourcePath, "w") as fw:
            fw.write(sHeaderText)
            fw.write("{:d}\n".format(self.n_internal_cells))
            fw.write("(\n")
            for val in self.reaction_rates:
                fw.write("{:e}\n".format(val))
            fw.write(")\n;\n\n")

            fw.write(boundaryText)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                fw.write("        type            fixedValue;\n")
                valString = "        value          uniform 0;\n"
                fw.write(valString)
                fw.write("    }\n")

            fw.write(closerText)

        if activation_dataset_error != "":
            zeroSourceErrorPath = zeroFolder / "SourceError"

            eHeaderText = """
        FoamFile
        {
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0";
            object      SourceError;
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        dimensions      [0 0 -1 0 0 0 0];

        internalField   nonuniform List<scalar>
        """
            with open(zeroSourceErrorPath, "w") as fw:
                fw.write(eHeaderText)
                fw.write("{:d}\n".format(self.n_internal_cells))
                fw.write("(\n")
                for val in self.sampled_rr_error:
                    fw.write("{:e}\n".format(val))
                fw.write(")\n;\n\n")

                fw.write(boundaryText)

                for face in self.patches.values():
                    fw.write("    " + face.face_id + "\n    {\n")
                    fw.write("        type            fixedValue;\n")
                    valString = "        value          uniform 0;\n"
                    fw.write(valString)
                    fw.write("    }\n")

                fw.write(closerText)

        return

    def generate_tr_source_file(self, time_treatment):
        """
        this file create the source file for the time residency fictiotious
        scalar in the zero folder.
        There are two modes:
        1) steady state: the source is 1 in every cell
        2) transient the source is 0 in every cell
        """

        zeroFolder = self.path / "0"
        zeroSourcePath = zeroFolder / "TrSource"

        sHeaderText = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      TrSource;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

"""
        intFieldText = """

internalField   uniform {};

        """
        boundaryText = """
boundaryField
{

"""
        closerText = """

}


// *********************************************************************** //
"""

        with open(zeroSourcePath, "w") as fw:
            fw.write(sHeaderText)

            if time_treatment == "steadystate":
                fw.write(intFieldText.format(1))
            elif time_treatment == "transient":
                fw.write(intFieldText.format(0))

            fw.write(boundaryText)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                fw.write("        type            fixedValue;\n")
                valString = "        value          uniform 0;\n"
                fw.write(valString)
                fw.write("    }\n")

            fw.write(closerText)

        return

    def assign_isotope_data(self):
        """
        Using the decay constant, understand if we are considering N-16, N-17, O-19
        or F-20 and assign the spectrum accordingly.

        Logic:
        - If "custom" is set in the isotope section of the transportProperties file,
        fetch data from the local database (legacy workflow). Optionally also load
        OpenMC photon decay distribution if OpenMC is available.
        - Otherwise, fetch photon decay distribution from the OpenMC chain XML file.
        (Requires OpenMC.)
        """
        import math

        # Decay constants (1/s)
        N16_DECAY = 0.09721559
        N17_DECAY = 0.1661825
        O19_DECAY = 0.02578672546
        F20_DECAY = 0.062093271

        try:
            import openmc  # type: ignore[import]
        except ImportError:
            openmc = None

        if self.isotope == "custom":
            # Identify isotope from decay constant
            if math.isclose(self.decay_constant, N16_DECAY, rel_tol=1e-3):
                self.isotope = "n16"
            elif math.isclose(self.decay_constant, N17_DECAY, rel_tol=1e-3):
                self.isotope = "n17"
            elif math.isclose(self.decay_constant, O19_DECAY, rel_tol=1e-3):
                self.isotope = "o19"
            elif math.isclose(self.decay_constant, F20_DECAY, rel_tol=1e-3):
                self.isotope = "f20"
            else:
                raise ValueError("ERROR, isotope not recognized")

            isotope_database = load_isotopes()
            iso_key = self.isotope  # already lowercase in this branch
            if iso_key not in isotope_database:
                raise ValueError("ERROR isotope not found in the database")

            isotope_data = isotope_database[iso_key]
            self.e_lines = isotope_data.e_lines
            self.p_lines = isotope_data.p_lines
            self.e_bins = isotope_data.e_bins
            self.p_bins = isotope_data.p_bins
            self.tot_p_emission = isotope_data.tot_p_emission
            self.particle_type = isotope_data.emitting_particle

            # Optional OpenMC distribution (photon) if available
            if openmc is not None:
                self.isotope_openmc_format = iso_key[:1].upper() + iso_key[1:].lower()
                self.decay_energy_distribution_ev = openmc.data.decay_photon_energy(
                    self.isotope_openmc_format
                )

        else:
            # This path currently supports photon-emitting isotopes via OpenMC only.
            if openmc is None:
                raise ImportError(
                    "openmc is required to fetch isotope data from chain xml file"
                )

            iso = (self.isotope or "").strip()
            if not iso:
                raise ValueError("ERROR, isotope is empty")

            self.isotope_openmc_format = iso[:1].upper() + iso[1:].lower()

            # Stored in electronvolt
            self.decay_energy_distribution_ev = openmc.data.decay_photon_energy(
                self.isotope_openmc_format
            )

            # OpenMC provides energies in eV -> convert to MeV
            self.e_lines = self.decay_energy_distribution_ev.x / 1e6

            # Convert from rate-like intensity to photons/decay (divide by decay constant, 1/s)
            self.p_lines = (
                self.decay_energy_distribution_ev.p
                / openmc.data.decay_constant(self.isotope_openmc_format)
            )

            self.particle_type = "photon"
            self.e_bins, self.p_bins = bins_from_lines(self.e_lines, self.p_lines)
            self.tot_p_emission = sum(self.p_lines)

    def read_t(self):
        """
        this function reads the T scalar
        """

        print("reading scalar values...")

        # common patterns
        # internalBlockPat = re.compile(
        #     r"internalField.*?\((.{1,}?)\)", re.MULTILINE | re.DOTALL
        # )
        #
        # tFile = self.path / str(self.last_time) / "T"
        #
        # try:
        #     inpFile = open(tFile, "r", encoding="utf8", errors="ignore")
        # except IOError:
        #     raise FileNotFoundError("couldn't open T file")
        # with inpFile:
        #     text = inpFile.read()
        #     numInternalBlocks = internalBlockPat.findall(text)
        #     internalScalar = numInternalBlocks[0].split("\n")[1:-1]
        #
        #     self.t_scalar = np.zeros(self.n_internal_cells)
        #     for i in range(self.n_internal_cells):
        #         self.t_scalar[i] = float(internalScalar[i])

        self.t_scalar = self.foamlib_object[-1]["T"].internal_field

        return

    def write_sampled_cartesian_source_vtk(self):
        """
        this function write the sampled source in a vtk file
        """

        out_vtk_file_path = self.results_folder / "cartesian_sampled_source.vtk"

        dims = (self.x_ints + 1, self.y_ints + 1, self.z_ints + 1)

        spacing = (
            (self.x_nodes[-1] - self.x_nodes[0]) / self.x_ints,
            (self.y_nodes[-1] - self.y_nodes[0]) / self.y_ints,
            (self.z_nodes[-1] - self.z_nodes[0]) / self.z_ints,
        )

        origin = (self.x_nodes[0], self.y_nodes[0], self.z_nodes[0])

        raw = np.asarray(
            [
                vox["normalized_emission_rate_per_voxel"]
                for vox in self.cartesian_voxel_list
            ],
            dtype=float,
        )
        reordered = raw.reshape(
            (self.x_ints, self.y_ints, self.z_ints), order="C"
        ).ravel(order="F")

        dataset_name_output = "emission_rate_per_voxel"

        write_cartesian_vtk(
            out_vtk_file_path, dims, spacing, origin, reordered, dataset_name_output
        )

        return

    def write_cdgs(self):
        """
        this function write the sample CDGS file
        """

        print("writing source model file in CDGS format ...")

        cdgsFile = self.results_folder / "cartesian_sampled_source.cdgs"

        with open(cdgsFile, "w") as fw:
            fw.write("num_meshes 1\n")
            fw.write("global_source {:e}\n".format(self.total_isotope_emission_rate))
            fw.write("mesh_id 1\n")
            fw.write("cdgs from vtk\n")
            fw.write("Cooling_time 0.0\n")
            fw.write("total_source {:e}\n".format(self.total_isotope_emission_rate))
            fw.write("energy_type {}\n".format("bins"))
            fw.write("energy_boundaries {:d}\n".format(len(self.e_bins)))
            # WRITE SPECTRUM BINS
            specString = formatValues(self.e_bins)
            fw.write(specString)

            fw.write("mesh_type rec\n")
            fw.write(
                "mesh_boundaries {:d} {:d} {:d}\n".format(
                    self.x_ints + 1, self.y_ints + 1, self.z_ints + 1
                )
            )
            fw.write("0.000000e+00  0.000000e+00  0.000000e+00\n")
            fw.write("1.000000e+00  0.000000e+00  0.000000e+00\n")
            fw.write("0.000000e+00  1.000000e+00  0.000000e+00\n")
            xString = formatValues(self.x_nodes)
            fw.write(xString)
            yString = formatValues(self.y_nodes)
            fw.write(yString)
            zString = formatValues(self.z_nodes)
            fw.write(zString)
            fw.write("source_data\n")

            voxelString1 = "{:d} {:.5e} {:.5e} 1\n"
            voxelString2 = "0 1.0 {:.5e}\n"
            specErrorString = formatValues([0] * (len(self.p_bins)))

            for vox in self.cartesian_voxel_list:
                if vox["normalized_emission_rate_per_voxel"] > 0:
                    fw.write(
                        voxelString1.format(
                            vox["id"],
                            vox["normalized_emission_rate_per_voxel"],
                            self.cartesian_voxel_volume_cm3,
                        )
                    )
                    fw.write(
                        voxelString2.format(vox["normalized_emission_rate_per_voxel"])
                    )
                    emittingSpectrum = [
                        val * vox["normalized_emission_rate_per_voxel"]
                        for val in self.p_bins
                    ]
                    spectrumString = formatValues(emittingSpectrum)
                    fw.write(spectrumString)
                    fw.write(specErrorString)

            fw.write("end_source_data\n")

        return

    def write_openmc_sm_source(self):
        """
        this function write the mesh-based radiation source file
        """

        print(
            "writing openmc source model file of the FLUNED results sampled over a cartesian grid  ..."
        )

        import openmc  # type: ignore[import]

        if openmc is None:
            print("ERROR openmc not available in the environment")
            print("skipping writing of the radiation source model.")
            return

        mesh_lower_left = [min(self.x_nodes), min(self.y_nodes), min(self.z_nodes)]

        mesh_upper_right = [max(self.x_nodes), max(self.y_nodes), max(self.z_nodes)]

        openmc_source_file_name = "settings_sm_fluned_source.xml"
        settings_source_file = self.results_folder / openmc_source_file_name

        source_mesh = openmc.RegularMesh()
        source_mesh.lower_left = mesh_lower_left
        source_mesh.upper_right = mesh_upper_right
        source_mesh.dimension = (self.x_ints, self.y_ints, self.z_ints)

        raw_strengths = np.asarray(
            [
                voxel["normalized_emission_rate_per_voxel"]
                for voxel in self.cartesian_voxel_list
            ],
            dtype=float,
        )
        # RegularMesh source strengths are indexed with x varying fastest.
        strengths = raw_strengths.reshape(
            (self.x_ints, self.y_ints, self.z_ints), order="C"
        ).ravel(order="F")

        source_mesh_space_dist = openmc.stats.MeshSpatial(
            mesh=source_mesh,
            strengths=strengths,
            volume_normalized=False,
        )

        rad_source = openmc.IndependentSource()
        rad_source.angle = openmc.stats.Isotropic()
        rad_source.energy = self.decay_energy_distribution_ev
        rad_source.space = source_mesh_space_dist
        rad_source.particle = self.particle_type
        rad_source.strength = self.total_isotope_emission_rate

        settings = openmc.Settings()
        settings.run_mode = "fixed source"
        settings.batches = 10
        settings.particles = int(1e5)
        settings.photon_transport = True
        settings.source_rejection_fraction = 1e-4
        settings.source = rad_source
        settings.export_to_xml(
            path=str(settings_source_file),
        )

        return settings_source_file

    def write_openmc_um_source(self):
        """
        This function write the unstructured mesh-based radiation source file
        based on the vtk file - generates a settings.xml with the mesh specification
        """
        print("writing openmc unstructured mesh source file ..")

        import openmc  # type: ignore[import]

        if openmc is None:
            print("ERROR openmc not available in the environment")
            print("skipping writing of the radiation source model.")
            return

        openmc.config["resolve_paths"] = False

        openmc_source_file_name = "settings_um_fluned_source.xml"
        self.settings_source_file = self.results_folder / openmc_source_file_name

        h5m_basename = "um_geometry.h5m"
        self.mesh_source_file = self.results_folder / h5m_basename

        vtk_tri_mesh_source_file = self.results_folder / "triangularized_t_scalar.vtk"

        # compute values and store in memory
        self.compute_triangularized_emission_rates(
            vtk_tri_mesh_source_file, scaling_factor=100.0, save_tri_mesh_vtk=True
        )

        # save a scaled up h5m file (m -> cm)
        generate_triangularized_h5m_um_mesh(self.vtk_file_path, self.mesh_source_file)

        #
        source_mesh = openmc.UnstructuredMesh(
            filename=str(h5m_basename), library="moab"
        )

        source_mesh_space_dist = openmc.stats.MeshSpatial(
            mesh=source_mesh,
            strengths=self.tri_mesh_emission_rates,
            volume_normalized=False,
        )

        rad_source = openmc.IndependentSource()
        rad_source.angle = openmc.stats.Isotropic()
        rad_source.energy = self.decay_energy_distribution_ev
        rad_source.space = source_mesh_space_dist
        rad_source.particle = self.particle_type
        rad_source.strength = self.total_isotope_emission_rate

        settings = openmc.Settings()
        settings.run_mode = "fixed source"
        settings.batches = 10
        settings.particles = int(1e5)
        settings.photon_transport = True
        settings.source_rejection_fraction = 1e-4
        settings.source = rad_source
        settings.export_to_xml(
            path=str(self.settings_source_file),
        )

        return self.settings_source_file

    # def write_openmc_um_source(self, mesh_id=100):
    #     """
    #     This function write the unstructured mesh-based radiation source file
    #     based on the vtk file
    #     """
    #
    #     print("writing openmc unstructured mesh source file ..")
    #
    #     openmc_source_file_name = "um_fluned_source.xml"
    #     openmc_source_file = self.results_folder / openmc_source_file_name
    #
    #     h5m_basename = "um_geometry.h5m"
    #     openmc_source_mesh_file = self.results_folder / h5m_basename
    #
    #     vtk_tri_mesh_source_file = self.results_folder / "triangularized_t_scalar.vtk"
    #
    #     openmc_source_import_commands = self.results_folder / "openmc_um_commands.txt"
    #
    #     self.compute_triangularized_emission_rates(
    #         vtk_tri_mesh_source_file, scaling_factor=100.0, save_tri_mesh_vtk=True
    #     )
    #
    #     # save a scaled up h5m file
    #     generate_triangularized_h5m_um_mesh(self.vtk_file_path, openmc_source_mesh_file)
    #
    #     with open(openmc_source_import_commands, "w") as f:
    #         f.write("from lxml import etree\n")
    #         f.write("parser = etree.XMLParser(huge_tree=True)\n")
    #         f.write(
    #             'source_root = etree.parse("{}", parser=parser).getroot()\n'.format(
    #                 openmc_source_file_name
    #             )
    #         )
    #         f.write("mesh_element = source_root.find('mesh')\n")
    #         f.write(
    #             "mesh_geo = openmc.UnstructuredMesh.from_xml_element(mesh_element)\n"
    #         )
    #         f.write("source_element = source_root.find('source')\n")
    #         f.write(
    #             "source = openmc.IndependentSource.from_xml_element(source_element, {100:mesh_geo})\n"
    #         )
    #
    #     # create root element
    #     root = et.Element("source")
    #
    #     # create sublement with the mesh source
    #     source_mesh = et.SubElement(
    #         root,
    #         "source",
    #         type="independent",
    #         particle=self.particle_type,
    #         strength=f"{self.total_isotope_emission_rate:.6e}",
    #     )
    #
    #     space = et.SubElement(
    #         source_mesh,
    #         "space",
    #         type="mesh",
    #         mesh_id=str(mesh_id),
    #         volume_normalized="False",
    #     )
    #     strengths = et.SubElement(space, "strengths")
    #     strengths.text = " ".join(
    #         map(
    #             str,
    #             [val for val in self.tri_mesh_emission_rates],
    #         )
    #     )  # adjust that the decay rate has been calculated after scaling the vtk file
    #
    #     # angle = et.SubElement(source_mesh, "angle", type="isotropic")
    #
    #     ebins_temp = [e * 1e6 for e in self.e_lines]  # convert from MeV to eV
    #     energy_parameters = [*ebins_temp, *self.p_lines]
    #     energy = et.SubElement(source_mesh, "energy", type="discrete")
    #     params = et.SubElement(energy, "parameters")
    #     params.text = " ".join(map(str, energy_parameters))
    #
    #     # create mesh element with id attribute
    #     mesh = et.SubElement(
    #         root,
    #         "mesh",
    #         id=str(mesh_id),
    #         name="source_mesh",
    #         type="unstructured",
    #         library="moab",
    #     )
    #
    #     # add child elements with text content
    #     filename = et.SubElement(mesh, "filename")
    #     filename.text = h5m_basename
    #
    #     # write to file with xml declaration
    #     tree = et.ElementTree(root)
    #     tree.write(
    #         openmc_source_file,
    #         encoding="utf-8",
    #         pretty_print=True,
    #         xml_declaration=True,
    #     )
    #
    #     return

    def compute_triangularized_emission_rates(
        self,
        tri_mesh_vtk_filepath,
        scaling_factor=100.00,
        cell_data_array="T",
        save_tri_mesh_vtk=True,
    ):
        """
        this function takes the vtk of the simulation, triangularizes it,
        extracts the T array, volumes and computes the emission rates
        (decay particle per second)
        """

        # a scaling factor is applied to generate the mesh for the source with the
        # units in cm for subsequent simulation

        generate_triangularized_scalar_mesh(
            self.vtk_file_path, tri_mesh_vtk_filepath, scaling_factor, cell_data_array
        )

        # extract the data from the triangularized mesh
        tri_mesh_volumes = get_vtk_volumes(tri_mesh_vtk_filepath)
        tri_mesh_isotopes_coonc = get_vtk_celldata_array(
            tri_mesh_vtk_filepath, cell_data_array
        )

        # the 1e6 factor is to take into account that the concentration data is taken by the scaled mesh
        # switch from emission per m3 to per cm3
        self.tri_mesh_emission_rates = [
            conc * vol * self.decay_constant * self.tot_p_emission / (scaling_factor**3)
            for conc, vol in zip(tri_mesh_isotopes_coonc, tri_mesh_volumes)
        ]

        tot_triangularized = sum(self.tri_mesh_emission_rates)

        if not math.isclose(
            tot_triangularized,
            self.total_isotope_emission_rate,
            rel_tol=1e-5,
            abs_tol=0,
        ):
            print("ERROR consistency issue after traingularizing the array")

            print(
                "total emission rate from tetras in the vtk file: ",
                sum(self.tri_mesh_emission_rates),
            )
            print(
                "total emission rate from the vtk file: ",
                self.total_isotope_emission_rate,
            )

            print()

            raise ValueError("ERROR consistency issue after traingularizing the array")

        if not save_tri_mesh_vtk:
            Path(tri_mesh_vtk_filepath).unlink()

        return

    def generate_system_files(self, time_treatment, fv_scheme):
        """
        this function creates the files needed in the system folder for an
        openFOAM simulation - in later development it will apply the case
        parameters
        """

        control_dict_text_transient = """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      controlDict;
    }

    application    FLUNED-solver;

    startFrom       startTime;

    startTime       0;

    stopAt          endTime;

    endTime         100;

    deltaT          0.1;

    writeControl    timeStep;

    writeInterval   100;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;

    adjustTimeStep  yes;

    maxCo           0.5;

    includeDecayScalar true;




    functions
    {

    """

        control_dict_text = """
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      controlDict;
    }

    application    FLUNED-solver;

    startFrom       startTime;

    startTime       0;

    stopAt          endTime;

    endTime         1000;

    deltaT          1;

    writeControl    timeStep;

    writeInterval   100;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;

    includeDecayScalar true;




    functions
    {

    """
        vol_calc_text = """
        volumeCalc
        {{
            type            writeCellVolumes;
            libs            ("libfieldFunctionObjects.so");
            select      all;

            writeFields     false;
            writeControl {};

        }}

    """

        vol_flow_text = """
        {}
        {{

            type            surfaceFieldValue;
            libs            ("libfieldFunctionObjects.so");
            patch   {};
            fields  (phi);

            operation sum;
            select  patch;
            name        $patch;

            writeFields     false;
            writeControl {};

        }}

    """
        vol_tx_flow_text = """
        {}
        {{

            type            surfaceFieldValue;
            libs            ("libfieldFunctionObjects.so");
            patch   {};
            fields  ({});
            weightField phi;

            operation sum;
            select  patch;
            name        $patch;

            writeFields     false;
            writeControl {};

        }}

    """

        vol_tx_sum_text = """
        {}
        {{

            type            volFieldValue;
            libs            ("libfieldFunctionObjects.so");
            fields  ({});

            operation       volAverage;
            select      all;

            writeFields     false;
            writeControl {};


        }}

    """
        system_folder = self.path / "system"
        control_dict_path = system_folder / "controlDict"

        with open(control_dict_path, "w", encoding="utf-8") as fw:
            if time_treatment == "steadystate":
                fw.write(control_dict_text)
                write_control = "outputTime"
            elif time_treatment == "transient":
                fw.write(control_dict_text_transient)
                write_control = "timeStep"
            else:
                raise ValueError("time_treatment argument not recognized")

            fw.write(vol_calc_text.format(write_control))

            fw.write(vol_tx_sum_text.format("volTSum", "T", write_control))
            fw.write(vol_tx_sum_text.format("volTaSum", "Ta", write_control))
            fw.write(vol_tx_sum_text.format("volTdSum", "Td", write_control))

            for face in self.patches.values():
                if face.face_type in ["inlet", "outlet"]:
                    fw.write(
                        vol_flow_text.format(
                            "volFlow-" + face.face_id, face.face_id, write_control
                        )
                    )

                    fw.write(
                        vol_tx_flow_text.format(
                            "volTFlow-" + face.face_id,
                            face.face_id,
                            "T",
                            write_control,
                        )
                    )
                    fw.write(
                        vol_tx_flow_text.format(
                            "volTrFlow-" + face.face_id,
                            face.face_id,
                            "Tr",
                            write_control,
                        )
                    )
                    fw.write(
                        vol_tx_flow_text.format(
                            "volTdFlow-" + face.face_id,
                            face.face_id,
                            "Td",
                            write_control,
                        )
                    )
                    fw.write(
                        vol_tx_flow_text.format(
                            "volTaFlow-" + face.face_id,
                            face.face_id,
                            "Ta",
                            write_control,
                        )
                    )

            fw.write("}")

        fv_scheme_stable = """
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      fvSchemes;
    }}

    ddtSchemes
    {{
        default         {};
    }}

    gradSchemes
    {{
        default         cellLimited Gauss linear 1;
    }}

    divSchemes
    {{
        default         none;
        div(phi,T)      Gauss  upwind;
        div(phi,Ta)      Gauss upwind;
        div(phi,Td)      Gauss upwind;
        div(phi,Tr)      Gauss upwind;
    }}

    laplacianSchemes
    {{
        default         none;
        laplacian(DT,T) Gauss linear corrected;
        laplacian(Dturbulent,T) Gauss linear corrected;
        laplacian(DT,Tr) Gauss linear corrected;
        laplacian(Dturbulent,Tr) Gauss linear corrected;
        laplacian(DT,Ta) Gauss linear corrected;
        laplacian(Dturbulent,Ta) Gauss linear corrected;
        laplacian(DT,Td) Gauss linear corrected;
        laplacian(Dturbulent,Td) Gauss linear corrected;
    }}

    interpolationSchemes
    {{
        default         linear;
    }}

    snGradSchemes
    {{
        default         limited 0.5;
    }}

    """

        fv_scheme_accurate = """
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      fvSchemes;
    }}

    ddtSchemes
    {{
        default         {};
    }}

    gradSchemes
    {{
        default         cellLimited Gauss linear 0.5;
    }}

    divSchemes
    {{
        default         none;
        div(phi,T)      Gauss linearUpwind default;
        div(phi,Ta)      Gauss linearUpwind default;
        div(phi,Td)      Gauss linearUpwind default;
        div(phi,Tr)      Gauss linearUpwind default;
    }}

    laplacianSchemes
    {{
        default         none;
        laplacian(DT,T) Gauss linear corrected;
        laplacian(Dturbulent,T) Gauss linear corrected;
        laplacian(DT,Tr) Gauss linear corrected;
        laplacian(Dturbulent,Tr) Gauss linear corrected;
        laplacian(DT,Ta) Gauss linear corrected;
        laplacian(Dturbulent,Ta) Gauss linear corrected;
        laplacian(DT,Td) Gauss linear corrected;
        laplacian(Dturbulent,Td) Gauss linear corrected;
    }}

    interpolationSchemes
    {{
        default         linear;
    }}

    snGradSchemes
    {{
        default         limited 1;
    }}

    """

        schemes_path = system_folder / "fvSchemes"

        with open(schemes_path, "w", encoding="utf-8") as fw:
            # select the fv scheme
            if fv_scheme == "accurate":
                fv_scheme_text = fv_scheme_accurate
            elif fv_scheme == "stable":
                fv_scheme_text = fv_scheme_stable
            else:
                raise ValueError("fv_scheme argument not recognized")

            # select the time treatment
            if time_treatment == "steadystate":
                fw.write(fv_scheme_text.format("steadyState"))
            elif time_treatment == "transient":
                fw.write(fv_scheme_text.format("Euler"))

        fv_solution_text = """
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      fvSolution;
    }}
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    solvers
    {{
        T
        {{
            solver          PBiCGStab;
            preconditioner  DILU;
            tolerance       1e-07;
            relTol          {};
        }}
        Tr
        {{
            solver          PBiCGStab;
            preconditioner  DILU;
            tolerance       1e-07;
            relTol          {};
        }}
        Ta
        {{
            solver          PBiCGStab;
            preconditioner  DILU;
            tolerance       1e-07;
            relTol          {};
        }}
        Td
        {{
            solver          PBiCGStab;
            preconditioner  DILU;
            tolerance       1e-07;
            relTol          {};
        }}
    }}

    SIMPLE
    {{
        nNonOrthogonalCorrectors 0;
        consistent      yes;

        residualControl
        {{

            T               1e-6;
            Tr              1e-6;
            Ta              1e-6;
            Td              1e-6;

        }}


    }}

    """
        solution_path = system_folder / "fvSolution"

        with open(solution_path, "w", encoding="utf-8") as fw:
            if time_treatment == "steadystate":
                fw.write(fv_solution_text.format(0.01, 0.01, 0.01, 0.01))
            elif time_treatment == "transient":
                fw.write(fv_solution_text.format(0, 0, 0, 0))

        parallel_dict_text = """
    FoamFile
    {
        format      ascii;
        class       dictionary;
        object      decomposeParDict;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    numberOfSubdomains 4;

    method          scotch;

    // *********************************************************************** //
    """
        parallel_dict_path = system_folder / "decomposeParDict"
        with open(parallel_dict_path, "w", encoding="utf-8") as fw:
            fw.write(parallel_dict_text)

        return

    def generate_constant_file(
        self, isotope, molecular_diffusion, decay_constant, schmidt_number
    ):
        """
        This function creates the files in the constant folder
        """

        transport_prop_text = """
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      transportProperties;
    }}
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    transportModel  Newtonian;

    isotope          {};

    DT             DT [0 2 -1 0 0 0 0] {};

    lambda         lambda [0 0 -1 0 0 0 0] {};

    Sct            Sct [ 0 0 0 0 0 0 0 ] {};


    // ************************************************************************ //
    """

        constant_folder = self.path / "constant"
        transport_prop_path = constant_folder / "transportProperties"

        with open(transport_prop_path, "w", encoding="utf-8") as fw:
            fw.write(
                transport_prop_text.format(
                    isotope,
                    molecular_diffusion,
                    decay_constant,
                    schmidt_number,
                )
            )

        return

    def write_empty_openmc_model(
        self,
        *,
        openmc_model_folder: Path | str = "openmc_empty_model",
        radsource_settings_file: Path | None = None,
        copy_radsource_mesh: bool = False,
        meshtally_voxel_size: float = 1.0,
        meshtally_max_voxel_n: float = 1e6,
        padding_factor: float = 1.5,
    ):
        """
        write the openmc geometry of the simulation circuit
        """

        import openmc  # type: ignore[import]

        if openmc is None:
            raise RuntimeError("openmc is required")

        mat_filename = "fluned_materials.xml"
        openmc_model_folder = self.results_folder / openmc_model_folder
        openmc_model_folder.mkdir(parents=True, exist_ok=True)

        mat_full_path = openmc_model_folder / mat_filename
        mats = openmc.Materials()

        # mats.cross_sections = str(parameters["xsection_file_path"])

        mats.export_to_xml(path=str(mat_full_path))

        geo_filename = "fluned_geometry.xml"
        geo_full_path = openmc_model_folder / geo_filename

        openmc_cells = []

        # converts vtk dimensions from m to cm and applies a padding factor
        scale_factor = 100.0

        # scaled bounds
        x0 = self.vtk_dimensions[0] * scale_factor
        x1 = self.vtk_dimensions[1] * scale_factor
        y0 = self.vtk_dimensions[2] * scale_factor
        y1 = self.vtk_dimensions[3] * scale_factor
        z0 = self.vtk_dimensions[4] * scale_factor
        z1 = self.vtk_dimensions[5] * scale_factor

        # centers
        cx = 0.5 * (x0 + x1)
        cy = 0.5 * (y0 + y1)
        cz = 0.5 * (z0 + z1)

        # half-sizes (extents) scaled by padding_factor
        hx = 0.5 * (x1 - x0) * padding_factor
        hy = 0.5 * (y1 - y0) * padding_factor
        hz = 0.5 * (z1 - z0) * padding_factor

        lower_left = [cx - hx, cy - hy, cz - hz]
        upper_right = [cx + hx, cy + hy, cz + hz]

        px0 = openmc.XPlane(x0=lower_left[0], boundary_type="vacuum")
        px1 = openmc.XPlane(x0=upper_right[0], boundary_type="vacuum")
        py0 = openmc.YPlane(y0=lower_left[1], boundary_type="vacuum")
        py1 = openmc.YPlane(y0=upper_right[1], boundary_type="vacuum")
        pz0 = openmc.ZPlane(z0=lower_left[2], boundary_type="vacuum")
        pz1 = openmc.ZPlane(z0=upper_right[2], boundary_type="vacuum")

        universe_box = +px0 & -px1 & +py0 & -py1 & +pz0 & -pz1

        universe_cell = openmc.Cell(
            name="outside",
            fill=None,
            region=universe_box,
        )

        openmc_cells.append(universe_cell)

        geo = openmc.Geometry(openmc_cells)
        geo.export_to_xml(path=str(geo_full_path))

        tallies_filename = "fluned_tallies.xml"
        tallies_full_path = openmc_model_folder / tallies_filename

        mesh = openmc.RegularMesh()

        mesh.lower_left = lower_left
        mesh.upper_right = upper_right

        mesh.dimension = mesh_ints_from_bounds(
            lower_left,
            upper_right,
            meshtally_voxel_size,
            int(meshtally_max_voxel_n),
        )

        # LARGE MESHTALLIES WILL QUICKLY CRASH THE WSL
        mesh_filter = openmc.MeshFilter(mesh)
        tally = openmc.Tally(name="flux_meshtally")
        tally.filters = [mesh_filter]
        tally.scores = ["flux"]
        tallies = openmc.Tallies([tally])

        tallies.export_to_xml(path=str(tallies_full_path))

        if radsource_settings_file is not None:
            shutil.copy(
                radsource_settings_file,
                openmc_model_folder,
            )
            if copy_radsource_mesh:
                shutil.copy(
                    self.mesh_source_file,
                    openmc_model_folder,
                )

            with warnings.catch_warnings():
                copied_settings_path = (
                    openmc_model_folder / Path(radsource_settings_file).name
                )
                warnings.simplefilter("ignore", openmc.IDWarning)
                settings = openmc.Settings().from_xml(copied_settings_path)

        else:
            # completely empty settings
            settings_filename = "settings.xml"
            settings_full_path = openmc_model_folder / settings_filename
            settings = openmc.Settings()
            settings.export_to_xml(path=str(settings_full_path))
        #
        model_filename = "model.xml"
        model_full_path = openmc_model_folder / model_filename

        model = openmc.Model(
            geometry=geo, materials=mats, settings=settings, tallies=tallies
        )
        model.export_to_model_xml(path=str(model_full_path))

        return
