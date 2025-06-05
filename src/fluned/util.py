import linecache
import tracemalloc
import os
import sys
import gzip
import vtk
from vtk.util import numpy_support as VN


def check_int(check_str):
    """
    Check if a string can be converted to an integer
    """
    try:
        int(check_str)
        return True
    except ValueError:
        return False


def mod(vec1):
    sum = 0
    for el in vec1:
        sum += el**2

    return pow(sum, 0.5)


def check_float(str):
    """
    this function checks if a string is a float
    """
    try:
        float(str)
        return True
    except ValueError:
        return False


def open_utf8_or_gzip(file_path):
    """
    this function tries to open a file with utf-8 encoding, if it fails it
    tries to open it with gzip. It returns the file object
    """

    try:
        with open(file_path, "rb") as f:
            magic = f.read(2)
    except OSError:
        print("Couldn't open Volume V file")
        sys.exit(1)

    if magic == b"\x1f\x8b":
        # The file is gzip-compressed
        try:
            with gzip.open(file_path, "rt", encoding="utf-8") as inpFile:
                data = inpFile.read()
        except Exception as e:
            print(f"Error reading gzip file: {e}")
            sys.exit(1)

    else:
        # The file is a regular text file
        try:
            with open(file_path, "r", encoding="utf-8") as inpFile:
                data = inpFile.read()
        except Exception as e:
            print(f"Error reading text file: {e}")
            sys.exit(1)

    return data


def sample_coordinates_vtk(vtk_file, dataset_name, coordinates):
    """
    this function reads the vtk and sample some datasets
    """

    check_file = os.path.isfile(vtk_file)
    if not check_file:
        print("ERROR activation file not found")
        sys.exit()

    # read the vtk file with an unstructured grid
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(vtk_file)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()

    # define probe

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(coordinates))

    # print ("parsed  0/{:d}".format(len(coordinates)))
    # step = pow(10,int(math.log10(len(coordinates)))-1)

    for i, val in enumerate(coordinates):
        points.SetPoint(i, val[0], val[1], val[2])
    #     if (i+1) % step == 0:
    #         print ("parsed  {:d}/{:d}".format(i+1,len(coordinates)))

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)

    # Perform the interpolation
    probe_filter = vtk.vtkProbeFilter()
    probe_filter.SetSourceData(data)
    probe_filter.SetInputData(polydata)
    probe_filter.Update()

    vtk_array = probe_filter.GetOutput().GetPointData().GetArray(dataset_name)

    reac_rates = VN.vtk_to_numpy(vtk_array)

    return reac_rates


def display_top(snapshot, key_type="lineno", limit=10):
    snapshot = snapshot.filter_traces(
        (
            tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
            tracemalloc.Filter(False, "<unknown>"),
        )
    )
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        print(
            "#%s: %s:%s: %.1f KiB"
            % (index, frame.filename, frame.lineno, stat.size / 1024)
        )
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print("    %s" % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))

    return

def formatValues(vector):
    maxLen = 70
    returnString = ""
    newLine = ""

    for item in vector:
        newNumber = "{:.7e}".format(item)
        if returnString == "" and newLine == "":
            newLine = newNumber
            continue

        if len(newLine + " " + newNumber) > maxLen:
            returnString += newLine + "\n"
            newLine = newNumber

        else:
            newLine = newLine + " " + newNumber

    returnString += newLine + "\n"

    return returnString
