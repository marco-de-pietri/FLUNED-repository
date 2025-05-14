# -*- coding: utf-8 -*-
import re
import sys
import os
import copy
import math
import numpy as np
import argparse
import linecache
import tracemalloc
import vtk
import pyvista as pv
from vtk.util import numpy_support as VN
import gzip
from water_isotopes.water_isotopes import get_isotope_data
from flundePost import flundePost

__version__ = "0.0.1"

def main():


    parser = argparse.ArgumentParser(description="flunedPost version " + __version__)
    parser.add_argument("-v","--vtk",  type=str,
            help="vtk file to sample", default = False)
    parser.add_argument("--cdgs", action ='store_true' ,
            help="print cdgs", default = False)
    parser.add_argument("--openmc", action ='store_true' ,
            help="generate openmc mesh-based source", default = False)
    parser.add_argument("--precision", type=float, default=0.01,
            help="sampling precision in m")
    parser.add_argument("--dataset", type=str, default='T',
            help="dataset in the vtk file to sample")
    parser.add_argument("--scaling", type=float, default=1,
            help="scaling factor for the total emission")
    args = parser.parse_args()




    simulationFolder = os.getcwd()
    simCase = flunedPostCase(simulationFolder)
    simCase.parseConstants()
    simCase.getIsotope()
    simCase.getTimeTreatment()
    simCase.get_last_time_step()
    simCase.getVTKPath(args.vtk)

    simCase.read_volumes()

    if args.source_cdgs or args.source_openmc:

        simCase.precision = args.precision
        simCase.dataset = args.dataset
        simCase.scaling = args.scaling
        simCase.getOriginalEmission()
        simCase.calculateSamplingCoordinates()
        simCase.sampleCoordinatesValues()

        if args.cdgs:
            simCase.writeCDGS()

        if  args.openmc:
            simCase.writeOpenMCSource()

if __name__ == '__main__':
    main()