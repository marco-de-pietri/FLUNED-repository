# -*- coding: utf-8 -*-
import argparse

from .. import __version__
from ..fluned_case_class import flunedCase
from ..input_handler import create_input_template, read_fluned_input_file


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
        raise ValueError("No input provided, printintg template and exiting")

    if args.input_template:
        print("printing template and exiting")
        create_input_template()
        return

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

        print("copying openfoam files from last CFD iteration ... ")
        fCase.generate_fluned_files()

        if args.launch_simulation:
            fCase.launch_solver()

        print("FINISHED!")

    return


if __name__ == "__main__":
    main()
