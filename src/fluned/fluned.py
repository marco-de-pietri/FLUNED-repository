# -*- coding: utf-8 -*-
import argparse
import copy
import os
import pickle
import re
import sys

import h5py

from . import __version__
from .fluned_case_class import flunedCase
from .openmc_chain_file_util import (
    filter_channels,
    map_targets_to_channels,
)


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

        # print(channels)

        new_case["openmc_depletion_channels"] = channels

        cases.append(new_case)

    return cases


def get_isotope_reactions_dicts(file_path):
    """
    this function  reads the depletion results and gets the dictionaries mapping the reactions and the parent isotopes
    """

    with h5py.File(file_path, "r") as f:
        nuc_idx = pickle.loads(f["index_nuc"][...].tobytes())
        rx_idx = pickle.loads(f["index_rx"][...].tobytes())

    return nuc_idx, rx_idx


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
TIME_TREATMENT  steadyState          #steadyState or transient supported
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


def read_fluned_input_file(path):
    """
    this function reads the user input file
    """

    case_path = re.compile(
        r"^\s*case .{1,}?(?=^\s*case|\Z)", re.MULTILINE | re.DOTALL | re.IGNORECASE
    )
    cases_vec = []
    parameters = [
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
        "div_scheme",
        "cfd_type",
        "fluent_fluid_region_name",
    ]
    try:
        with open(path, "r", encoding="utf8", errors="ignore") as fin:
            text_block = fin.read()
    except IOError:
        print("couldn't open file")
        text_block = ""

    cases_blocks = case_path.findall(text_block)

    for case in cases_blocks:
        parameters_dict = {}
        case_lines = case.splitlines()
        for line in case_lines:
            if len(line.strip()) == 0:
                continue
            if '"' in line:
                args = line.strip().split('"')
                key = args[0].strip().lower()
                if key in parameters:
                    parameters_dict[key] = args[1]
            else:
                args = line.strip().split()
                if args[0].lower() in parameters:
                    parameters_dict[args[0].lower()] = args[1]

        if (
            "simulation_type" in parameters_dict
            and parameters_dict["simulation_type"].lower() == "openmc-multi"
        ):
            openmc_cases_dicts = generate_openmc_simulation_parameters(parameters_dict)
            cases_vec.extend(openmc_cases_dicts)
        else:
            cases_vec.append(parameters_dict)

    return cases_vec


def main():
    parser = argparse.ArgumentParser(
        description="FLUNED case generator version " + __version__
    )
    parser.add_argument("-i", "--input", type=str, help="input")
    parser.add_argument(
        "-l",
        "--launch_simulation",
        action="store_true",
        help="launch simulation",
        default=False,
    )
    parser.add_argument(
        "-t",
        "--input_template",
        action="store_true",
        help="create an input template",
        default=False,
    )
    args = parser.parse_args()

    if not args.input and not args.input_template:
        print("WARNING no input provided")
        print("printing template and exiting")
        create_input_template()
        sys.exit()

    if args.input_template:
        print("printing template and exiting")
        create_input_template()
        sys.exit()

    input_cases = read_fluned_input_file(args.input)

    # define a vector of FLUNED cases

    for case in input_cases:
        print("creating FLUNED case...")

        fCase = flunedCase()
        fCase.generate_case(case)
        fCase.initialize_cfd_class()

        if fCase.cfd_type == "fluent-h5-multi":
            print("parsing h5 files ... ")
            fCase.parse_fluent_simulation()
            fCase.create_case_folders_fluent()
            fCase.convert_fluent_to_openfoam()

        if fCase.cfd_type == "fluent-multi":
            raise NotImplementedError(
                "Parsing binary cas/dat files not implemented yet"
            )

        print("copying openfoam files from last CFD iteration ... ")
        fCase.generate_fluned_files()

        if args.launch_simulation:
            fCase.launch_solver()

        print("FINISHED!")

    return


if __name__ == "__main__":
    main()
