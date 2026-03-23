# -*- coding: utf-8 -*-
import argparse
import tracemalloc
from pathlib import Path

from .. import __version__
from ..fluned_case_class import flunedCase
from ..utils import display_top


def main():
    parser = argparse.ArgumentParser(description="flunedPost version " + __version__)
    parser.add_argument(
        "-d", "--debug", action="store_true", help="debug mode", default=False
    )
    parser.add_argument(
        "-c", "--check", action="store_true", help="check mode", default=False
    )
    parser.add_argument(
        "-s", "--cdgs", action="store_true", help="print cdgs", default=False
    )
    parser.add_argument(
        "--openmc-sm",
        action="store_true",
        help="generate a openmc regular mesh source format",
        default=False,
    )
    parser.add_argument(
        "--openmc-um",
        action="store_true",
        help="generate a openmc unstructured mesh source format",
        default=False,
    )
    parser.add_argument(
        "--precision", type=float, default=1, help="sampling precision in cm"
    )
    parser.add_argument(
        "--dataset", type=str, default="T", help="dataset in the vtk file to sample"
    )
    args = parser.parse_args()

    if args.debug:
        tracemalloc.start()

    simulation_folder = Path.cwd()
    fluned_case = flunedCase(fluned_path=simulation_folder)
    fluned_case.parse_fluned_simulation()

    if args.check:
        fluned_case.launch_check_mode()

    # write source files
    if args.cdgs or args.openmc_sm or args.openmc_um:
        fluned_case.source_sampling_resolution_cm = args.precision
        fluned_case.source_sampling_dataset = args.dataset

        if args.cdgs or args.openmc_sm:
            fluned_case.fluned_simulation.run_cartesian_sampling(
                fluned_case.source_sampling_resolution_cm,
                fluned_case.source_sampling_dataset,
            )

        if args.cdgs:
            fluned_case.fluned_simulation.write_cdgs()

        if args.openmc_sm:
            settings_file = fluned_case.fluned_simulation.write_openmc_sm_source()
            fluned_case.fluned_simulation.write_empty_openmc_model(
                openmc_model_folder="structured_mesh_source",
                radsource_settings_file=settings_file,
                copy_radsource_mesh=False,
            )

        if args.openmc_um:
            settings_file = fluned_case.fluned_simulation.write_openmc_um_source()
            fluned_case.fluned_simulation.write_empty_openmc_model(
                openmc_model_folder="unstructured_mesh_source",
                radsource_settings_file=settings_file,
                copy_radsource_mesh=True,
            )

    fluned_case.write_results(args)

    print("Finished!")

    if args.debug:
        print(tracemalloc.get_traced_memory())
        snapshot = tracemalloc.take_snapshot()
        display_top(snapshot)

    return


if __name__ == "__main__":
    main()
