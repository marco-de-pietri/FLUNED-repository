# -*- coding: utf-8 -*-
import re
import sys
import os
import argparse
from .fluned_case_class import flunedCase

__version__ = "0.1.0"


def create_input_template():
    """
    this function generates an input template file
    """

    input_name = "inputTemplate"
    current_folder = os.getcwd()

    input_path = os.path.join(current_folder, input_name)

    template_text = """CASE  FLUNED_01_DEFAULT_N16
TIME_TREATMENT  steadyState #steadyState or transient supported
#ACTIVATION_FILE  {}
#ACTIVATION_DATASET    "Value - Total"
#ACTIVATION_CONST 1e16
#ACTIVATION_NORMALIZATION 0 #leave zero if no normalization is required
#DECAY_CONSTANT 0.0257867           #O19 decay const
#DECAY_CONSTANT 0.1661825           #N17 decay const
DECAY_CONSTANT 0.09721559           #N16 decay const
INLET_CONC 1e10
MOLECULAR_DIFFUSION 2e-09
SCHMIDT_NUMBER   0.7
CFD_PATH      "{}"
CFD_TYPE    OpenFoam       # OpenFoam, fluent-h5-multi types supported
#FLUENT_FLUID_REGION_NAME     region_name
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
        fin = open(path, "r", encoding="utf8", errors="ignore")
    except IOError:
        print("couldn't open file")
    with fin:
        text_block = fin.read()

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

        if fCase.cfd_type == "fluent-h5-multi":
            print("parsing h5 files ... ")
            fCase.getH5files()
            fCase.create_case_folders_fluent()
            fCase.getFluidCells_h5()
            fCase.getFluidFaces_h5()
            fCase.read_density_h5()
            fCase.defineWalls_multi_h5()
            fCase.writeOwner_multi_h5()
            fCase.writeNeighbour_multi_h5()
            fCase.write_cell_zones()
            fCase.writeBoundary_multi_h5()
            fCase.writeFaces_multi_h5()
            fCase.writeNodes_multi_h5()
            fCase.writePhi_multi_h5()
            fCase.writeSpeed_multi_h5()
            fCase.writeNut_multi_h5()

        if fCase.cfd_type == "fluent-multi":
            print("parsing binary cas/dat files not implemented yet ... ")
            sys.exit()

        print("copying last CFD iteration files ... ")
        fCase.create_case_folder()
        fCase.initialize_cfd_class()
        fCase.copy_last_phi()
        fCase.copy_last_u()
        fCase.getNumCells()
        fCase.copyLastNut()
        fCase.copyPolyMesh()
        fCase.reconstructFaces()
        fCase.generate_system_file()
        fCase.launch_vol_func_object()
        fCase.generate_constant_file()
        fCase.generate_zero_t()
        fCase.generate_zero_ta()
        fCase.generate_zero_td()
        fCase.generate_zero_tr()
        fCase.generateSourceFile()
        fCase.generateTrSourceFile()

        if args.launch_simulation:
            fCase.launch_solver()

        print("FINISHED!")

    return


if __name__ == "__main__":
    main()
