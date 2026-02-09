# -*- coding: utf-8 -*-
import argparse
import os

from . import __version__
from .fluned_case_multi_isotopes_class import flunedCaseMultiIsotopes


def main():
    parser = argparse.ArgumentParser(
        description="fluned-post-multi-isotopes " + __version__
    )
    parser.add_argument(
        "-f",
        "--folder",
        type=str,
        default=os.getcwd(),
        help="folder containing the fluned simulations",
    )
    parser.add_argument(
        "--dataset", type=str, default="T", help="dataset in the vtk file to sample"
    )
    args = parser.parse_args()

    multi_isotopes_case = flunedCaseMultiIsotopes(args.folder)

    multi_isotopes_case.post_process_cases()

    multi_isotopes_case.generate_openmc_um_source()

    print("Finished")

    return


if __name__ == "__main__":
    main()
