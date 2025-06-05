# -*- coding: utf-8 -*-
import os
import argparse
import tracemalloc
from .util import display_top

from .fluned_case_class import flunedCase

__version__ = "0.1.0"


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
        "--openmc",
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
        "--precision", type=float, default=0.01, help="sampling precision in m"
    )
    parser.add_argument(
        "--dataset", type=str, default="T", help="dataset in the vtk file to sample"
    )
    parser.add_argument(
        "--scaling", type=float, default=1, help="scaling factor for the total emission"
    )
    args = parser.parse_args()

    if args.debug:
        tracemalloc.start()

    simulation_folder = os.getcwd()
    # simCase = flunedPostCase(simulationFolder)
    simCase = flunedCase(fluned_path=simulation_folder)
    simCase.parse_constant()
    simCase.get_isotope()
    simCase.get_time_treatment()
    simCase.get_last_time_step()

    if args.check:
        simCase.launch_grad_func_object()
        simCase.read_velocities()
        simCase.read_grad_t()

    simCase.read_volumes()
    simCase.read_t()
    simCase.readPostProcess_flows()
    simCase.readPostProcess_Tflows()
    simCase.readPostProcess_Trflows()
    simCase.generate_vtk()
    simCase.create_results_folder()

    # write source files
    if args.cdgs or args.openmc or args.openmc_um:
        simCase.precision = args.precision
        simCase.dataset = args.dataset
        simCase.scaling = args.scaling
        simCase.get_vtk_path()
        simCase.get_total_decays()
        simCase.write_external_stl()

        if args.cdgs or args.openmc:
            simCase.calculateSamplingCoordinates()
            simCase.sampleCoordinatesValues()
            simCase.write_sampled_source_vtk()

        if args.cdgs:
            simCase.writeCDGS()

        elif args.openmc:
            simCase.writeOpenMCSource()

        elif args.openmc_um:
            simCase.write_openmc_um_source()

    simCase.write_summary(args)

    print("Finished!")

    if args.debug:
        print(tracemalloc.get_traced_memory())
        snapshot = tracemalloc.take_snapshot()
        display_top(snapshot)

    return


if __name__ == "__main__":
    main()
