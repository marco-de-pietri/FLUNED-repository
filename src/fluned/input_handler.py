import copy
import re
from pathlib import Path

import numpy as np
import pyvista as pv

from .multi_isotope_utils import (
    filter_channels,
    get_isotope_reactions_dicts,
    map_targets_to_channels,
)

# these are parameters that can be specified for each isotope in multi-isotope simulations.
_INPUT_PARAMETERS_ISOTOPE = [
    "inlet_conc",
    "schmidt_number",
    "molecular_diffusion",
]

_INPUT_PARAMETERS_REAL_LIST = [
    "activation_rotation_degs",
    "activation_translation_m",
]
_INPUT_PARAMETERS = [
    "case",
    "time_treatment",
    "simulation_type",
    "activation_constant",
    "activation_dataset",
    "activation_dataset_error",
    "activation_rotation_center_mode",
    "activation_rotation_euler_order",
    "fv_scheme",
    "isotope",
    "activation_file",
    "activation_normalization",
    "inlet_conc",
    "decay_constant",
    "cfd_path",
    "molecular_diffusion",
    "schmidt_number",
    "cfd_type",
    "fluent_fluid_region_name",
]


def fluned_defaults():
    # Canonical defaults
    return {
        "simulation_type": "single-isotope",
        "fv_scheme": "stable",
        "molecular_diffusion": 2e-09,
        "decay_constant": 0.0,
        "isotope": "custom",
        "activation_file": "",
        "activation_dataset": "",
        "activation_dataset_error": "",
        "activation_normalization": 0.0,
        "activation_constant": 0.0,
        "activation_rotation_center_mode": "origin",
        "activation_rotation_euler_order": "xyz",
        "activation_rotation_degs": np.array([0.0, 0.0, 0.0]),
        "activation_translation_m": np.array([0.0, 0.0, 0.0]),
    }


def get_openmc_reaction_channels(chain_file_path, isotopes_index, reactions_index):
    """
    This function takes the openmc chain file and returns the channel information that
    leads to the formation of concerning radioisotopes
    """

    all_targets = map_targets_to_channels(
        chain_file_path, isotopes_index, reactions_index
    )

    emitting_targets = filter_channels(chain_file_path, all_targets)

    return emitting_targets


def get_activation_dataset_name(
    vtk_file: str | Path, test_dataset_name: str | None = None
) -> str | None:
    """
    Open a VTK file with PyVista and inspect cell data arrays.
    If exactly one cell-data array exists, return (name, array).
    Otherwise return None.
    """
    mesh = pv.read(str(vtk_file))

    keys = list(mesh.cell_data.keys())
    if test_dataset_name is not None:
        if test_dataset_name in keys:
            return test_dataset_name
        else:
            return None

    if len(keys) != 1:
        return None

    name = keys[0]
    return name


def get_single_vtk_file(folder: str | Path) -> Path:
    """
    Return the Path of the single .vtk file in the folder.
    If there are zero or more than one .vtk files, return Path().
    """
    folder_path = Path(folder)
    vtk_files = list(folder_path.glob("*.vtk"))

    return vtk_files[0].resolve() if len(vtk_files) == 1 else Path()


def create_input_template():
    """
    this function generates an input template file
    """

    input_name = "input_template"
    current_folder = Path.cwd()

    input_path = current_folder / input_name
    comment_string = "# "
    dataset_name = "Value - Total"

    vtk_file = get_single_vtk_file(current_folder)

    if not vtk_file.is_dir():
        activation_comment = ""
        vtk_path = vtk_file.as_posix()
        vtk_dataset_name = get_activation_dataset_name(vtk_path)
        if vtk_dataset_name:
            dataset_name = vtk_dataset_name
    else:
        activation_comment = comment_string
        vtk_path = current_folder.as_posix() + "/path/to/activation_file.vtk"

    template_text = f"""CASE  FLUNED_01_DEFAULT
TIME_TREATMENT  steadyState          # steadyState or transient supported
SIMULATION_TYPE single-isotope       # single-isotope, openmc-multi supported
{activation_comment}ACTIVATION_FILE  {vtk_path}
{activation_comment}ACTIVATION_DATASET    "{dataset_name}"
# ACTIVATION_CONSTANT 1e16
# ACTIVATION_NORMALIZATION 0         # Leave zero if no normalization is required
INLET_CONC 1e10
ISOTOPE    <name>
SCHMIDT_NUMBER   0.7
CFD_PATH      "{current_folder.as_posix()}"
CFD_TYPE    OpenFoam                 # OpenFoam, fluent-h5-multi types supported
# FLUENT_FLUID_REGION_NAME     region_name
"""

    input_path.write_text(template_text, encoding="utf-8")

    return


def generate_openmc_simulation_parameters(args_dict):
    """
    this function parse the depletion file and the chains.xml and return a dictionary
    with the simulation parameters for the openmc cases, one for each isotope
    """

    cases = []

    import openmc  # type: ignore[import]

    depletion_isotopes, depletion_reactions = get_isotope_reactions_dicts(
        args_dict["activation_file"]
    )

    target_isotopes_channels = get_openmc_reaction_channels(
        openmc.config["chain_file"], depletion_isotopes, depletion_reactions
    )

    for target, channels in target_isotopes_channels.items():
        new_case = copy.deepcopy(args_dict)

        isotope = target.lower()

        new_case["isotope"] = isotope
        new_case["case"] = new_case["case"] + "_" + target.upper()
        new_case["openmc_depletion_channels"] = channels

        # if it was defined it overrides the generic argument
        for parameter in _INPUT_PARAMETERS_ISOTOPE:
            isotope_specific_parameter = f"{parameter}_{isotope}"
            if isotope_specific_parameter in args_dict:
                new_case[parameter] = args_dict[isotope_specific_parameter]

        cases.append(new_case)

    return cases


def read_fluned_input_file(path):
    """
    this function reads the user input file
    """

    case_path = re.compile(
        r"^\s*case .{1,}?(?=^\s*case|\Z)", re.MULTILINE | re.DOTALL | re.IGNORECASE
    )
    cases_vec = []
    try:
        with open(path, "r", encoding="utf8", errors="ignore") as fin:
            text_block = fin.read()
    except IOError:
        print("couldn't open file")
        text_block = ""

    cases_blocks = case_path.findall(text_block)
    defaults = fluned_defaults()

    for case in cases_blocks:
        parameters_dict = {}
        case_lines = case.splitlines()
        for line in case_lines:
            if len(line.strip()) == 0:
                continue
            if '"' in line:
                args = line.strip().split('"')
                key = args[0].strip().lower()
                if key in _INPUT_PARAMETERS:
                    parameters_dict[key] = args[1]
            else:
                args = line.strip().split()
                if args[0].lower() in _INPUT_PARAMETERS:
                    parameters_dict[args[0].lower()] = args[1]
                elif args[0].lower().rsplit("_", 1)[0] in _INPUT_PARAMETERS_ISOTOPE:
                    parameters_dict[args[0].lower()] = args[1]
                if args[0].lower() in _INPUT_PARAMETERS_REAL_LIST:
                    parameters_dict[args[0].lower()] = np.array(
                        [float(x) for x in args[1:]]
                    )

        # Apply defaults (user values override defaults)
        parameters_dict = {**defaults, **parameters_dict}

        if parameters_dict.get("simulation_type", "").lower() == "openmc-multi":
            openmc_cases_dicts = generate_openmc_simulation_parameters(parameters_dict)
            cases_vec.extend(openmc_cases_dicts)
        else:
            cases_vec.append(parameters_dict)

    return cases_vec
