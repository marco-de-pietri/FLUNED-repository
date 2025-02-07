# -*- coding: utf-8 -*-
"""
this script calculates the activation in a given circuit
"""
import argparse
import os
import sys

from .utils import create_input_template, read_input_file

from .fluned_sl_case import FlunedSlCase

__version__ = "0.1.0"

def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(description="FLUNED_SL - version: " + __version__)
    parser.add_argument("-i", "--input", type=str, help="input")
    parser.add_argument(
        "-t",
        "--input_template",
        action="store_true",
        help="create an input template",
        default=False,
    )
    parser.add_argument(
        "-s",
        "--sample_rrmesh",
        action="store_true",
        help="sample the rr mesh before running the simulation",
        default=False,
    )
    args = parser.parse_args()

    if args.input_template:
        print("printing template and exiting")
        create_input_template()
        sys.exit()

    if not args.input:
        print("WARNING no input provided")
        print("printing template and exiting")
        create_input_template()
        sys.exit()

    if args.sample_rrmesh:
        print("sampling the rr mesh and exiting")
        user_values = read_input_file(args.input)
        case_path = os.getcwd()
        case = FlunedSlCase(user_values, case_path)



    user_values = read_input_file(args.input)

    case_path = os.getcwd()

    case = FlunedSlCase(user_values, case_path)

    case.solve()

    # case.print_all_nodes()

    # list_out = case.get_outflow_node_list()
    # tot_in = case.get_total_inflow()
    # tot_out = case.get_total_outflow()
    # print (tot_in)
    # print (list_out)
    # print (tot_out)

    # test traverseNodes

    # test1 = ("01_CIRCUIT1", 369)

    # traversed_nodes_list = case.traverse_nodes(test1)
    # traversed_nodes_list = list(reversed(traversed_nodes_list))

    # for node in traversed_nodes_list:
    #    print (node)

    return 0


if __name__ == "__main__":
    main()
