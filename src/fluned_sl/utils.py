"""
this file contains some utility functions
"""
import os
import re
import sys
from distutils.util import strtobool
import pyvista as pv
import numpy as np
import vtk
from concurrent.futures import ThreadPoolExecutor


def interp_series(x_vals, y_vals, x_new):
    """
    this function returns a y value for a given x value, by interpolation of a
    series of x and y values
    """
    y_new = 0

    if len(x_vals) != len(y_vals):
        print("ERROR: x and y vectors must have the same length")
        sys.exit()

    if x_new < x_vals[0] or x_new > x_vals[-1]:
        y_new = 0
        return y_new

    for i, x_val in enumerate(x_vals):
        if x_new == x_val:
            y_new = y_vals[i]
            return y_new

        if x_new > x_vals[i] and x_new < x_vals[i+1]:
            y_new = ( y_vals[i] + ( (x_new - x_vals[i] ) *
                                    (y_vals[i + 1] - y_vals[i])  /
                                    (x_vals[i + 1] - x_vals[i]))
                     )
            return y_new

    return y_new

def bin_search(arr, target):
    """
    Performs a binary search to find the closest value to a target value in a sorted list.
    Returns the closest value.
    """
    left = 0
    right = len(arr) - 1

    while left <= right:
        mid = (left + right) // 2

        if arr[mid] == target:
            return mid

        elif arr[mid] < target:
            left = mid + 1

        else:
            right = mid - 1

    # Handle the case where the target value is not in the list
    if right < 0 or left >= len(arr):
        return -1

    # Compare the target value to the values of the closest element and its two neighboring elements
    #if target - arr[right] < arr[left] - target:
    #    return arr[right]
    #else:
    #    return arr[left]

    return min(left,right)

def linear_interp(x_list,y_list, closest_left_idx, new_x, outside = 'max'):
    """
    Performs a linear interpolation of a series knowing the closest left
    index to the value to be interpolated. outside can be given to decide the
    treatment if the value is outside the list
    """

    new_y = 0


    if closest_left_idx == -1:
        if outside == 'max':
            new_y = max(y_list)
        if outside == 'zero':
            new_y = 0
        if outside == 'last':
            new_y = y_list[-1]
    else:
        x_1 = x_list[closest_left_idx]
        x_2 = x_list[closest_left_idx + 1]
        y_1 = y_list[closest_left_idx]
        y_2 = y_list[closest_left_idx + 1]
        new_y = (y_1 +
                     (new_x - x_1)*
                     (y_2 - y_1)/
                     (x_2 - x_1))

    return new_y


def unpack_mcnp(string):
    """
    this function creates a list of float values from the compact string
    notation used in MCNP
    """

    char_vector = string.split()

    float_vector = []

    pat_number = re.compile(r"^(\d+)i\Z", re.IGNORECASE)
    j = 0


    for i,element in enumerate(char_vector):

        read_error = False

        if is_float(element):
            coordinate = float(element)

            if j != 0:
                if coordinate > float_vector[j-1]:
                    float_vector.append(coordinate)
                    j += 1
                else:
                    print ("ERROR: coordinates are not sequential!")
                    read_error = True
                    sys.exit()
            else:
                float_vector.append(coordinate)
                j += 1

        else:
            interp_number = pat_number.findall(element)

            if len(interp_number) == 1:
                if i==0 or i == (len(char_vector)-1):
                    print("ERROR: vector starts or end with interpolation!")
                    read_error = True
                    sys.exit()
                else:
                    begin_point = float(char_vector[i-1])
                    end_point = float(char_vector[i+1])
                    if end_point < begin_point:
                        print ("ERROR: coordinates are not sequential!")
                        sys.exit()
                    else:

                        interp = int(interp_number[0])
                        step = (end_point - begin_point) / (interp + 1)
                        for k in range(interp):
                            new_value = begin_point + step*(k+1)
                            float_vector.append(new_value)
                            j += 1


            else:
                print ("ERROR: wrong format")
                sys.exit()


        if read_error:
            sys.exit()

    return float_vector

def is_float(num):
    """
    helper function
    """
    try:
        float(num)
        return True
    except ValueError:
        return False

def split_text_file_lines(path):
    """
    this function can be applied to open and split the lines of the
    following files: inlets.dat, links.dat, nodes.dat

    Args:
        path (string): path of the text file to operate

    Returns:
        list of string : cleaned up list with the string lines of the file
    """

    try:
        fin = open(path, "r", encoding="utf8", errors="ignore")
    except IOError:
        print(f"couldn't open {path}")
        sys.exit()
    with fin:
        text_block = fin.read()

    text_lines = text_block.splitlines()

    text_lines_split = []

    for line in text_lines:
        text_lines_split.append(re.split("[ ,]+", line.strip(", ")))

    return text_lines_split


def get_circuit_files(path):
    """
    this function returns the paths of the circuit files

    Args:
        path (string): path  of a folder

    Returns:
        dict: dictionary with the paths of the circuit files
    """

    inlets_pattern = re.compile(r".*inlets.dat\Z", re.IGNORECASE)
    links_pattern  = re.compile(r".*links.dat\Z", re.IGNORECASE)
    nodes_pattern  = re.compile(r".*nodes.dat\Z", re.IGNORECASE)
    mcnp_pattern   = re.compile(r".*mcnp.dat\Z", re.IGNORECASE)
    rr_mesh_pattern = re.compile(r"^.*rrmesh\.(?:vtk|vtr)$", re.IGNORECASE)

    file_list = os.listdir(path)

    file_paths = {}

    for file in file_list:
        inlets_files = inlets_pattern.findall(file)
        links_files  = links_pattern.findall(file)
        nodes_files  = nodes_pattern.findall(file)
        mcnp_files   = mcnp_pattern.findall(file)
        rrmesh_files = rr_mesh_pattern.findall(file)

        if len(inlets_files) == 1:
            file_paths["inletsPath"] = os.path.join(
                os.path.abspath(path), inlets_files[0]
            )

        if len(links_files) == 1:
            file_paths["linksPath"] = os.path.join(
                os.path.abspath(path), links_files[0]
            )

        if len(nodes_files) == 1:
            file_paths["nodesPath"] = os.path.join(
                os.path.abspath(path), nodes_files[0]
            )

        if len(mcnp_files) == 1:
            file_paths["mcnpPath"] = os.path.join(
                os.path.abspath(path), mcnp_files[0]
            )
        if len(rrmesh_files) == 1:
            file_paths["rrmeshPath"] = os.path.join(
                os.path.abspath(path), rrmesh_files[0]
            )

    # print (filePaths)

    return file_paths

def search_circuits(case_directory):
    """
    this function search for folders that contain the three
    files that describe a circuit
    *nodes.dat,*links.dat,*inlets.dat,*mcnp.dat
    """

    dir_list = os.listdir(case_directory)

    temp_list = []
    for item in dir_list:
        item_folder = os.path.join(case_directory, item)
        if os.path.isdir(item_folder):
            temp_list.append((item, item_folder))
    folder_list = []
    for item in temp_list:
        circuit_files = get_circuit_files(item[1])
        found_files = len(circuit_files)
        if 0 < found_files < 3:
            print("WARNING incomplete circuit description")
            print(item)
            sys.exit()
        if found_files in [3,4,5]:
            folder_list.append(item)

    if len(folder_list) == 0:
        print("ERROR no circuit found")
        sys.exit()

    return folder_list

def create_input_template():
    """
    this function creates a template of the input file
    """

    input_name = "input_template"
    current_folder = os.getcwd()

    input_path = os.path.join(current_folder, input_name)

    template_text = """
STEADY_STATE       TRUE        # false for transient studies
ISOTOPE             N16        # isotope to be studied
N_SAMPLING          1e6        # number of histories
FLUID_CARRIER     WATER        # only water is supported for now
NUMERICAL_METHOD DETERMINISTIC # monte_carlo or deterministic
    """

    with open(input_path, "w", encoding="utf8") as file_w:
        file_w.write(template_text)

    return 0


def read_input_file(path):
    """
    this function reads the user input file

    Args:
        path (string): path of the input file

    Returns:
        dict: dictionary with the user input values
    """
    parameters = ["steady_state",
                  "isotope",
                  "n_sampling_min",
                  "n_sampling_max",
                  "rrmesh_sampling_error_max",
                  "rrmesh_sampling_cm",
                  "rrmesh_sampling_field_name",
                  "rrmesh_sampling_errfield_name",
                  "rrmesh_scaling_factor",
                  "fluid_carrier",
                  "pipe_time_default",
                  "tank_time_default",
                  "mc_error_max",
                  "t_delta_sec",
                  "t_begin_sec",
                  "t_end_sec",
                  "t_bins",
                  "n_parallel_processes",
                  "numerical_method",
                  ]

    try:
        fin = open(path, "r", encoding="utf8", errors="ignore")
    except IOError:
        print("couldn't open file")
        sys.exit()
    with fin:
        text_block = fin.read()

    parameters_dic = {}
    case_lines = text_block.splitlines()
    for line in case_lines:
        if len(line.strip()) == 0:
            continue
        if '"' in line:
            args = line.strip().split('"')
            key = args[0].strip().lower()
            if key in parameters:
                if is_float(args[1]):
                    parameters_dic[key] = float(args[1])
                else:
                    parameters_dic[key] = args[1]
        else:
            args = line.strip().split()
            if args[0].lower() in parameters:
                if is_float(args[1]):
                    parameters_dic[args[0].lower()] = float(args[1])
                elif args[1].lower() in ["true", "false"]:
                    parameters_dic[args[0].lower()] = bool(strtobool(args[1].lower()))
                else:
                    parameters_dic[args[0].lower()] = args[1].lower()


    return parameters_dic

def extend_pnt(point,axis, distance):
    """
    this function returns a point by extending a point, along a vector,
    by a given distance
    """

    newpoint = [point[0] + axis[0]*distance,
                point[1] + axis[1]*distance,
                point[2] + axis[2]*distance]

    return newpoint



def vtk_sampling_0(file_path, coords, parameters):
    """
    Samples cell-based values from a VTK file at the given coordinates.
    For each coordinate, the cell that contains the point is identified using
    PyVista's find_containing_cell method. If the point is within a cell,
    the function retrieves the cell's data for a specified field and error.
    If the error is less than or equal to the maximum allowed error, the field's
    value is returned; otherwise, np.nan is returned. If no cell contains the
    coordinate, np.nan is returned.

    Args:
        file_path (str): Path to the VTK file.
        coords (list or array-like): List of (x, y, z) coordinates.
        parameters (dict): Dictionary with keys:
            - "rrmesh_sampling_error_max": maximum allowable error.
            - "rrmesh_sampling_field_name": name of the cell data field to sample.
            - "rrmesh_sampling_errfield_name": name of the error field in cell data.

    Returns:
        list: A list of sampled field values, or np.nan when sampling fails or
              the error exceeds the allowed maximum.
    """
    max_sample_error = parameters["rrmesh_sampling_error_max"]
    sampling_field = parameters["rrmesh_sampling_field_name"]
    error_field = parameters["rrmesh_sampling_errfield_name"]

    mesh = pv.read(file_path)
    results = []

    for coord in coords:
        cell_id = mesh.find_containing_cell(coord)
        if cell_id > 0:
            # Retrieve the cell's field and error values using the cell_id.
            field_value = mesh.cell_data[sampling_field][cell_id]
            error_value = mesh.cell_data[error_field][cell_id]

            print ("debugging field_value",field_value)
            print ("debugging error_value",error_value)

            if error_value <= max_sample_error:
                results.append(field_value)
            else:
                results.append(np.nan)

    print ("debugging results",results)

    mean = np.nanmean(results)

    return mean


def _find_cell_value(coord, locator, mesh, sampling_field, error_field, max_sample_error, tol, max_cell_size):
    pcoords = [0.0, 0.0, 0.0]
    weights = [0.0] * max_cell_size
    generic_cell = vtk.vtkGenericCell()
    cell_id = locator.FindCell(coord, tol, generic_cell, pcoords, weights)
    if cell_id < 0:
        return np.nan
    # Retrieve the cell's data values.
    field_value = mesh.cell_data[sampling_field][cell_id]
    error_value = mesh.cell_data[error_field][cell_id]
    return field_value if error_value <= max_sample_error else np.nan

def vtk_sampling(file_path, coords, parameters, use_multithreading=True):
    """
    Samples cell-based values from a VTK file at given coordinates using a faster approach.
    This version builds a vtkStaticCellLocator once and (optionally) uses multithreading to
    process many coordinate queries in parallel.

    Args:
        file_path (str): Path to the VTK file.
        coords (list/array-like): List of (x, y, z) coordinates.
        parameters (dict): Contains:
            - "rrmesh_sampling_error_max": maximum allowed error.
            - "rrmesh_sampling_field_name": name of the field in cell data.
            - "rrmesh_sampling_errfield_name": name of the error field in cell data.
        use_multithreading (bool): If True, uses ThreadPoolExecutor to process points in parallel.

    Returns:
        list: Sampled field values (or np.nan if not found or error exceeds limit).
    """
    max_sample_error = parameters["rrmesh_sampling_error_max"]
    sampling_field   = parameters["rrmesh_sampling_field_name"]
    error_field      = parameters["rrmesh_sampling_errfield_name"]

    # Read the VTK file.
    mesh = pv.read(file_path)

    # Build a fast cell locator using vtkStaticCellLocator.
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(mesh)
    locator.BuildLocator()

    tol = 1e-6
    max_cell_size = int(mesh.GetMaxCellSize()) if hasattr(mesh, "GetMaxCellSize") else 100

    if use_multithreading:
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(_find_cell_value, coord, locator, mesh,
                                sampling_field, error_field,
                                max_sample_error, tol, max_cell_size)
                for coord in coords
            ]
            results = [future.result() for future in futures]
    else:
        results = [
            _find_cell_value(coord, locator, mesh,
                             sampling_field, error_field,
                             max_sample_error, tol, max_cell_size)
            for coord in coords
        ]


    if len(results) == 0 or np.all(np.isnan(results)):
        mean = 0.0
    else:
        mean = np.nanmean(results)
    return mean


