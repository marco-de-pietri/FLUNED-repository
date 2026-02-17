import copy
import os
import re

from .activation_file_utils import get_isotope_reactions_dicts
from .chain_file_utils import (
    filter_channels,
    map_targets_to_channels,
)

# these are parameters that can be specified for each isotope in multi-isotope simulations.
_INPUT_PARAMETERS_ISOTOPE = [
    "inlet_conc",
    "schmidt_number",
    "molecular_diffusion",
]
_INPUT_PARAMETERS = [
    "case",
    "time_treatment",
    "simulation_type",
    "activation_const",
    "activation_dataset",
    "activation_dataset_error",
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
    # Canonical defaults for *one* case
    return {
        "simulation_type": "single-isotope",
        "fv_scheme": "stable",
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


def create_input_template():
    """
    this function generates an input template file
    """

    input_name = "input_template"
    current_folder = os.getcwd()

    input_path = os.path.join(current_folder, input_name)

    template_text = """CASE  FLUNED_01_DEFAULT_N16
TIME_TREATMENT  steadyState          # steadyState or transient supported
SIMULATION_TYPE single-isotope       # single-isotope, openmc-multi supported
# ACTIVATION_FILE  {}
# ACTIVATION_DATASET    "Value - Total"
# ACTIVATION_CONST 1e16
# ACTIVATION_NORMALIZATION 0         # Leave zero if no normalization is required
# DECAY_CONSTANT 0.0257867           # O19 decay const
# DECAY_CONSTANT 0.1661825           # N17 decay const
DECAY_CONSTANT 0.09721559            # N16 decay const
INLET_CONC 1e10
MOLECULAR_DIFFUSION 2e-09
SCHMIDT_NUMBER   0.7
CFD_PATH      "{}"
CFD_TYPE    OpenFoam                 # OpenFoam, fluent-h5-multi types supported
# FLUENT_FLUID_REGION_NAME     region_name
"""

    with open(input_path, "w", encoding="utf-8") as fw:
        fw.write(template_text.format(current_folder, current_folder))

    return 0


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

        # Apply defaults (user values override defaults)
        parameters_dict = {**defaults, **parameters_dict}

        if parameters_dict.get("simulation_type", "").lower() == "openmc-multi":
            openmc_cases_dicts = generate_openmc_simulation_parameters(parameters_dict)
            cases_vec.extend(openmc_cases_dicts)
        else:
            cases_vec.append(parameters_dict)

    return cases_vec
