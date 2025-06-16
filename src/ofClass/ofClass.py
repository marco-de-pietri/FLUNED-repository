"""
class for the OF simulations, this can be used to parse and generate files
"""

import copy
import math
import os
import sys
import re
import pathlib

import numpy as np
import pyvista as pv
from fluned.util import open_utf8_or_gzip
from .fluned_tool_launchers import launch_centroid_func_object
from .fluned_tool_launchers import generate_vtk
from .fluned_vtk_utils import (
    get_vtk_dimensions,
    sample_vtk,
    generate_external_stl,
    write_cartesian_vtk,
)

from isotopes.isotopes import load_isotopes


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


def check_float(s):
    """
    check if a string can be converted to float
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def get_post_process_list(path, dir_prefix, face_name, col_name):
    """
    this function gets the post process flow from the simulation folder
    """

    dir_name = dir_prefix + face_name
    flow_folders = [itm for itm in os.listdir(path) if dir_name in itm]

    if len(flow_folders) != 1:
        raise ValueError("Error with the number of flow folders")

    flow_files = get_post_files(os.path.join(path, flow_folders[0]))

    area_m2 = extract_area_field(flow_files[0])

    time_lists = post_file_array(flow_files, "Time")
    flow_lists = post_file_array(flow_files, col_name)

    time_list_sorted, flow_list_sorted = merge_continue_runs(
        time_lists,
        flow_lists,
    )

    return flow_list_sorted, time_list_sorted, area_m2


def extract_area_field(file_path):
    """
    This function reads a postprocess file and
    if the area value is present it returns it
    """

    area = 0

    try:
        postFile = open(file_path, "r", encoding="utf8", errors="ignore")
    except IOError:
        raise FileNotFoundError("couldn't open postprocess file")
    with postFile:
        lines = postFile.readlines()
        for line in lines:
            line = line.replace("#", "")
            wrds = line.split()
            if "area" in wrds[0].lower():
                area = float(wrds[-1])
                break

    return area


def merge_continue_runs(time_lists, data_lists):
    """
    this function takes time series segments and join them in a single one -
    the overlapping sections are removed
    """

    if len(time_lists) == 0 or len(data_lists) == 0:
        raise ValueError("Error with the number of post process files")

    if len(time_lists) != len(data_lists):
        print("ERROR mismatch in length of data series")
        print("number of time lists: ", len(time_lists))
        print("number of data lists: ", len(data_lists))
        raise ValueError("Error with the number of post process files")

    tot_len_time_lists = sum([len(x) for x in time_lists])
    tot_len_data_lists = sum([len(x) for x in data_lists])

    if tot_len_time_lists != tot_len_data_lists:
        print("ERROR mismatch in length of data series")
        print("number of time points: ", tot_len_time_lists)
        print("number of data points: ", tot_len_data_lists)
        raise ValueError("Error with the number of post process data points")

    temp_list = list(zip(time_lists, data_lists))

    # time_lists_sorted = sorted(time_lists,
    #                           key=lambda x:x[0])

    sorted_temp_list = sorted(temp_list, key=lambda x: x[0][0])

    time_lists_sorted, data_lists_sorted = zip(*sorted_temp_list)
    # data_lists_sorted = sorted(data_lists,
    #                           key=lambda x:x[0])

    time_series = time_lists_sorted[0]
    data_series = data_lists_sorted[0]

    if len(time_lists) > 1:
        for i, data_list in enumerate(time_lists_sorted):
            time_series_temp = copy.deepcopy(time_series)
            if i == 0:
                continue
            appended = False
            for j, val in enumerate(time_series_temp):
                if val >= time_lists_sorted[i][0]:
                    time_series = time_series[0:j]
                    time_series.extend(data_list)
                    data_series = data_series[0:j]
                    data_series.extend(data_lists_sorted[i])
                    appended = True
                    break

            if not appended:
                time_series.extend(data_list)
                data_series.extend(data_lists_sorted[i])

    return time_series, data_series


def post_file_array(path_list, name):
    """
    this function extracts a generic array contained in the post
    processing file
    """

    array_list = []

    for f_path in path_list:
        array = []
        index = -1

        try:
            post_file = open(f_path, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open postprocess file")
        with post_file:
            lines = post_file.readlines()
            for line in lines:
                line = line.replace("#", "")

                wrds = line.split()

                if len(wrds) == 0:
                    continue

                if "time" in wrds[0].lower():
                    index = wrds.index(name)
                    continue

                if index != -1:
                    array.append(float(wrds[index]))

            array_list.append(array)

    return array_list


def get_post_files(path):
    """
    this function crawl the folder to reach the data in the post process file
    """

    file_paths = []
    folder_itms = os.listdir(path)
    folder_itms = sorted(folder_itms, reverse=False, key=lambda x: float(x))

    for fld in folder_itms:
        path_1 = os.path.join(path, fld)
        if os.path.isdir(path_1):
            file_itms = os.listdir(path_1)
            complete_path = os.path.join(path_1, file_itms[0])
            file_paths.append(complete_path)

    return file_paths


class PatchClass:
    """
    this class represents a patch of the mesh of an openfoam simulation
    """

    def __init__(self, patch_dict, path):
        """
        Constructs all the necessary attributes for the PatchClass object.

        Parameters:
        -----------
        patch_dict : dict
            A dictionary with the patch information.
        path : str
            The complete path to the simulation folder.
        """

        self.simulation_path = path
        self.post_process_path = os.path.join(self.simulation_path, "postProcessing")
        self.face_id = patch_dict["face_id"]
        self.face_elements_n = patch_dict["face_elements_n"]
        self.face_first_element = patch_dict["face_first_element"]
        self.face_type = patch_dict["type"]
        self.face_sum_phis = patch_dict["sum_phis"]
        self.post_process_flow = []
        self.post_process_time = []
        self.post_process_t_flow = []
        self.post_process_ta_flow = []
        self.post_process_td_flow = []
        self.post_process_tr_flow = []
        self.t_conc_atoms_m3 = []
        self.ta_conc_atoms_m3 = []
        self.td_conc_atoms_m3 = []
        self.tr_conc = []

        self.area_m2 = 0

    def post_process_face(self):
        """
        this function reads all the post process files and calculates the
        different scalar flows across the face
        """

        if self.face_type == "wall":
            return None

        self.post_process_flow, self.post_process_time, self.area_m2 = (
            get_post_process_list(
                self.post_process_path, "volFlow-", self.face_id, "sum(phi)"
            )
        )
        self.post_process_t_flow, _, _ = get_post_process_list(
            self.post_process_path, "volTFlow-", self.face_id, "sum(T)"
        )
        self.post_process_ta_flow, _, _ = get_post_process_list(
            self.post_process_path, "volTaFlow-", self.face_id, "sum(Ta)"
        )
        self.post_process_td_flow, _, _ = get_post_process_list(
            self.post_process_path, "volTdFlow-", self.face_id, "sum(Td)"
        )
        self.post_process_tr_flow, _, _ = get_post_process_list(
            self.post_process_path, "volTrFlow-", self.face_id, "sum(Tr)"
        )

        self.t_conc_atoms_m3 = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_t_flow, self.post_process_flow)
        ]
        self.ta_conc_atoms_m3 = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_ta_flow, self.post_process_flow)
        ]
        self.td_conc_atoms_m3 = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_td_flow, self.post_process_flow)
        ]
        self.tr_conc = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_tr_flow, self.post_process_flow)
        ]

        # check the sign of the flow
        if self.post_process_flow[-1] > 0 and self.face_type != "outlet":
            raise ValueError("Error with the flow sign")
        elif self.post_process_flow[-1] < 0 and self.face_type != "inlet":
            raise ValueError("Error with the flow sign")
        elif self.post_process_flow[-1] == 0 and self.face_type != "wall":
            raise ValueError("Error with the flow sign")
        if self.t_conc_atoms_m3[-1] < 0:
            raise ValueError("Error with the concentration sign")

        # print (self.post_process_t_flow)
        # print (self.t_conc_atoms_m3)

        return None


class SimulationOF:
    """
    A class to represent a simulation using OpenFOAM.
    some functions are present to parse and generate file

    Attributes:
    -----------
    path : str
        The complete path to the simulation folder.
    """

    def __init__(self, path: str):
        """
        Constructs all the necessary attributes for the SimulationOF object.

        Parameters:
        -----------
        path : str
            The complete path to the simulation folder.
        """
        self.path = path
        self.case = os.path.split(self.path)[-1]
        self.reduction_rate = []
        self.normalized_average_td = []
        self.inlet_td_conc_atoms_m3 = []
        self.outlet_rr_conc_atoms_m3 = []
        self.average_rr_conc_atoms_m3 = []
        self.post_process_td_average = []
        self.average_ta = []
        self.volume_m3 = 0
        self.vtk_file_path = ""
        self.scaled_vtk_file_path = ""
        self.vtk_dimensions = []

    def create_case_folder(self):
        """
        This function creates an OpenFoam Case - it calls some
        function for the definition of the specific files
        """

        # case Folder
        case_folder = self.path
        if not os.path.exists(case_folder):
            os.mkdir(case_folder)

        # T=0 folder
        zero_folder = os.path.join(case_folder, "0")
        if not os.path.exists(zero_folder):
            os.mkdir(zero_folder)

        # constant folder
        const_folder = os.path.join(case_folder, "constant")
        if not os.path.exists(const_folder):
            os.mkdir(const_folder)

        # polyMesh folder
        poly_folder = os.path.join(case_folder, "constant", "polyMesh")
        if not os.path.exists(poly_folder):
            os.mkdir(poly_folder)

        # system folder
        sys_folder = os.path.join(case_folder, "system")
        if not os.path.exists(sys_folder):
            os.mkdir(sys_folder)

        # case foam file
        case_file = os.path.join(case_folder, "case.foam")
        pathlib.Path(case_file).touch()

        return

    def post_process_openfoam_simulation(self):
        """
        this function is called when we process a finished openfoam simulation
        """

        self.post_process_path = os.path.join(self.path, "postProcessing")
        self.last_time = self.get_last_time()
        self.volumetric_flag = self.get_volumetric_flag()

        if not self.volumetric_flag:
            self.density_kg_m3 = self.get_density()
        else:
            self.density_kg_m3 = 1000

        self.patches = self.get_patches()
        self.n_internal_cells = self.get_number_internal_cells()

    def get_number_internal_cells(self):
        """
        this function create an attribute to the class that specifies the
        number of internal cells. It does so by reading the U file
        """

        velocityFile = os.path.join(self.path, str(self.last_time), "U")

        internal_block_pat = re.compile(
            r"internalField.*?(\d+).*?\(", re.MULTILINE | re.DOTALL
        )

        try:
            with open(velocityFile, "r", encoding="utf8", errors="ignore") as inp_file:
                text = inp_file.read()
        except (IOError, OSError) as e:
            raise RuntimeError(f"Unable to open velocity file '{velocityFile}': {e}")

        cell_number = internal_block_pat.findall(text)[0]

        return int(cell_number)

    def post_process_fluned_simulation(self):
        """
        this function is called when we process a finished FLUNED simulation
        """
        for face in self.patches.values():
            face.post_process_face()

        self.post_process_td_average, _, _ = get_post_process_list(
            self.post_process_path, "volTdSum", "", "volAverage(Td)"
        )
        self.average_ta, _, _ = get_post_process_list(
            self.post_process_path, "volTaSum", "", "volAverage(Ta)"
        )

        self.total_inlet_t_atoms = self.get_total_inlet_t_atoms()
        self.total_outlet_t_atoms = self.get_total_outlet_t_atoms()

        self.inlet_td_conc_atoms_m3 = self.get_inlet_td_conc_atoms_m3()
        self.outlet_t_conc_atoms_m3 = self.get_outlet_t_conc_atoms_m3()
        # print ("self.inlet_td_conc_atoms_m3")
        # print (f"{self.inlet_td_conc_atoms_m3[-1]:.2e}")

        self.reduction_rate = self.get_reduction_rate()
        # print ("self.reduction_rate")
        # print (self.reduction_rate)

        self.normalized_average_td = self.get_normalized_average_td()
        # print ("self.normalized_average_td")
        # print (self.normalized_average_td)

        self.outlet_rr_conc_atoms_m3 = self.get_outlet_rr_conc_atoms_m3()
        # print ("self.outlet_rr_conc_atoms_m3")
        # print (f"{self.outlet_rr_conc_atoms_m3[-1]:.2e}")

        # print ("self.average_ta")
        # print (f"{self.average_ta[-1]:.2e}")

        self.parse_constants_file()
        self.assign_isotope()
        self.get_time_treatment()
        self.read_volumes()
        self.read_t()

        self.create_results_folder()

        generate_vtk(self.path)
        self.vtk_file_folder = os.path.join(self.path, "VTK")
        self.get_vtk_file()
        self.vtk_dimensions, self.volume_m3 = get_vtk_dimensions(self.vtk_file_path)

        generate_external_stl(self.vtk_file_path, self.results_folder)

        self.isotope_amounts = [vol * t for vol, t in zip(self.volumes, self.t_scalar)]
        self.total_isotope_amount = sum(self.isotope_amounts)
        self.total_average_isotope_concentration = (
            self.total_isotope_amount / self.volume_m3
        )
        self.total_isotope_activity = self.decay_constant * self.total_isotope_amount
        self.total_isotope_emission_rate = (
            self.branching_ratio * self.total_isotope_activity
        )

        return None

    def get_vtk_file(self, custom_vtk=False):
        """
        this function looks for the vtk file in the VTK folder generated by the
        vtk openfoam utility and store it for generation of results

        custom values can be provided for special workflow where the vtk
        simulation file is rototranslated before generating the radiation source
        """

        vtkFiles = []

        self.vtk_path = ""

        if custom_vtk:
            filename = os.path.splitext(custom_vtk)[0]
        else:
            filename = ""

        pat_string = r"{}\.vtk\Z".format(filename)

        vtkFilePat = re.compile(pat_string, re.IGNORECASE)

        for filename in os.listdir(self.vtk_file_folder):
            matchVTKfiles = vtkFilePat.findall(filename)
            if len(matchVTKfiles) == 1:
                vtkFiles.append(filename)

        if len(vtkFiles) == 1:
            self.vtk_file_path = os.path.join(self.path, "VTK", vtkFiles[0])
        else:
            print("ERROR zero or more than one vtk files")
            print("found in the VTK folder")
            print("found vtk files: ", vtkFiles)
            raise ValueError()

        return

    def get_patches(self):
        """
        this function creates a dictionary with the patches objects
        """

        patches = {}

        patches_list = self.parse_boundary_phi_files()

        for patch in patches_list:
            patches[patch["face_id"]] = PatchClass(patch, self.path)

        return patches

    def get_reduction_rate(self):
        """
        this function calculates the reduction rate of the isotope in the cfd
        it sums all the outlet Td sum values - therefore multiple inlets or
        outlets are not supported yet
        """

        for patch in self.patches.values():
            if patch.face_type == "inlet":
                atom_in = patch.post_process_td_flow
            if patch.face_type == "outlet":
                atom_out = patch.post_process_td_flow

        red_ratio = [abs(x / y) if y != 0 else 0 for x, y in zip(atom_out, atom_in)]

        return red_ratio

    def get_outlet_rr_conc_atoms_m3(self):
        """
        this function returns the concentration of the outlet patch for the
        Ta field, meaning for the part generated only by the irradiation and
        not from the inlet flow
        at the moment it supports only one outlet patch
        """

        for patch in self.patches.values():
            if patch.face_type == "outlet":
                conc_out = patch.ta_conc_atoms_m3
                break

        return conc_out

    def get_inlet_td_conc_atoms_m3(self):
        """
        this function returns the concentration of the inlet patch for the
        Td field, meaning for the part due to the inlet flow

        at the moment it supports only one outlet patch
        """

        for patch in self.patches.values():
            if patch.face_type == "inlet":
                conc_out = patch.td_conc_atoms_m3
                break

        return conc_out

    def get_outlet_t_conc_atoms_m3(self):
        """
        this function returns the concentration of the outlet patch for the
        T field

        at the moment it supports only one outlet patch
        """

        for patch in self.patches.values():
            if patch.face_type == "outlet":
                conc_out = patch.t_conc_atoms_m3
                break

        return conc_out

    def get_total_inlet_t_atoms(self):
        """
        this function returns the total amount of atoms
        enetering the mesh through the inlets
        """

        atoms_in = 0

        for patch in self.patches.values():
            if patch.face_type == "inlet":
                atoms_in += patch.post_process_t_flow[-1]

        return atoms_in

    def get_total_outlet_t_atoms(self):
        """
        this function returns the total amount of atoms
        exiting the mesh through the inlets
        """

        atoms_out = 0

        for patch in self.patches.values():
            if patch.face_type == "outlet":
                atoms_out += patch.post_process_t_flow[-1]

        return atoms_out

    def get_normalized_average_td(self):
        """
        this function returns the average concentration of the Td field due
        only to the decay of the inlet flow
        """

        norm_td = [
            x / y if y != 0 else 0
            for x, y in zip(self.post_process_td_average, self.inlet_td_conc_atoms_m3)
        ]

        return norm_td

    def parse_boundary_phi_files(self):
        """
        this function examines phi files and returns a list of dictionary with
        the face_name and type
        at the moment it is not optimized for speed
        """

        # print ("parsing faces...")

        faces = []

        poly_mesh_folder = os.path.join(self.path, "constant", "polyMesh")
        boundary_file = os.path.join(poly_mesh_folder, "boundary")

        try:
            inp_file = open(boundary_file, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open boundary file")
        with inp_file:
            text = inp_file.read()

            face_def_pat = re.compile(
                r"\d+[\n\r\s]+?\(.*?[\n\r\s]+?\)", re.MULTILINE | re.DOTALL
            )
            face_number_pat = re.compile(
                r"(\d+)[\n\r\s]+?\(.*?\)", re.MULTILINE | re.DOTALL
            )
            boundary_pat = re.compile(
                r"[^\s]+[\n\r\s]+?\{.*?\}", re.MULTILINE | re.DOTALL
            )
            boundary_name_pat = re.compile(
                r"([^\s]+)[\n\r\s]+?\{.*?\}", re.MULTILINE | re.DOTALL
            )
            face_n_pat = re.compile(r"nFaces.*?(\d+)")
            first_face_pat = re.compile(r"startFace.*?(\d+)")
            def_block = face_def_pat.findall(text)[0]
            face_number = int(face_number_pat.findall(text)[0])
            boundary_blocks = boundary_pat.findall(def_block)

            for block in boundary_blocks:
                b_dict = {}
                b_dict["face_id"] = boundary_name_pat.findall(block)[0]
                b_dict["face_elements_n"] = int(face_n_pat.findall(block)[0])
                b_dict["face_first_element"] = int(first_face_pat.findall(block)[0])
                faces.append(b_dict)

            if len(faces) != face_number:
                raise ValueError("Error with the number of faces")

        phi_file_path = os.path.join(self.path, str(self.last_time), "phi")
        try:
            inp_file = open(phi_file_path, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open phi file")
        with inp_file:
            face_phi_pat = re.compile(r"\((.{1,}?)\)", re.MULTILINE | re.DOTALL)
            text = inp_file.read()
            wall_face_pat = re.compile(r"value\s+uniform\s+0")

            # print ("face phi")
            for face in faces:
                face_block_pat = re.compile(
                    face["face_id"] + r"[\n\r\s]+?\{.*?\}", re.MULTILINE | re.DOTALL
                )
                face_block = face_block_pat.findall(text)[0]
                # print (faceBloc)
                face_phis = face_phi_pat.findall(face_block)
                # print (facePhis)
                wall_confirm = wall_face_pat.findall(face_block)
                # print (wallConfirm)
                if len(face_phis) == 0 and len(wall_confirm) != 0:
                    face["type"] = "wall"
                    face["sum_phis"] = 0
                else:
                    phi_list = face_phi_pat.findall(face_block)[0].strip().split("\n")
                    phi_list = np.array([float(val) for val in phi_list])

                    if all(i >= 0 for i in phi_list):
                        face["type"] = "outlet"
                        face["sum_phis"] = sum(phi_list)
                    elif all(i <= 0 for i in phi_list):
                        face["type"] = "inlet"
                        face["sum_phis"] = sum(phi_list)
                    else:
                        raise ValueError("Error, phis with mixed sign in boundary")

        # for ls in faces:
        #    print (ls)

        return faces

    def get_last_time(self):
        """
        Returns the last time folder in the simulation folder.
        """

        time_folders = self.get_time_folders()

        return max(time_folders)

    def get_time_folders(self):
        """
        Returns a list of the time folders in the simulation folder.

        Returns:
        --------
        list
            A list of the time folders in the simulation folder.
        """
        time_folders = []
        # List all entries in the directory
        entries = os.listdir(self.path)

        for name in entries:
            if not os.path.isdir(os.path.join(self.path, name)):
                continue
            try:
                float(name)
                time_folders.append(name)
            except ValueError:
                continue
        return time_folders

    def get_volumetric_flag(self):
        """
        this function check the phi values to check if the phis are in the units
        of volumetric or mass flux (m3/s or kg/s)
        """

        dimension_vec = self.get_phi_dimensions()

        if dimension_vec == [1, 0, -1, 0, 0, 0, 0]:
            volumetric_flag = False
        elif dimension_vec == [0, 3, -1, 0, 0, 0, 0]:
            volumetric_flag = True
        else:
            raise ValueError(
                "The dimensions of phi are not recognized: ", dimension_vec
            )

        return volumetric_flag

    def get_density(self):
        """
        this function gets the density value from the
        constant/transportProperties file
        """

        file_path = self.path + "/constant/transportProperties"
        dict_path = "rhoRef"

        value = self.query_of_single_value(file_path, dict_path)

        return float(value)

    def get_phi_dimensions(self):
        """
        this function gets the dimensions of the phi values from the last
        time folder
        """

        file_path = self.path + "/" + str(self.last_time) + "/phi"

        dimension_vec = self.query_dimensions(file_path)

        return dimension_vec

    def query_dimensions(self, file_path):
        """
        this function queries the phi file to get the dimensions of the phi
        values
        """

        dict_path = "dimensions"

        dim_vector = []

        for tok_dict in self.token_classifier(file_path):
            if tok_dict["path"] == dict_path and tok_dict["type"] == "VECTOR_DATA":
                dim_vector.append(int(tok_dict["token"]))
            if len(dim_vector) == 7:
                break

        return dim_vector

    def query_of_single_value(self, file_path, dict_path, position=-2):
        """
        this function returns one value from one key of a OF dictionary

        Args:
            file_path (string): OF file path
            dict_path (string): dictionary path
            position (int, optional): token number of those with the same dict
                                      path. Defaults to -1.
        """
        field_vector = []

        found = False

        for tok_dict in self.token_classifier(file_path):
            if tok_dict["path"] == dict_path:
                field_vector.append(tok_dict["token"])
                found = True

            if found and tok_dict["path"] != dict_path:
                break

        # print (field_vector)

        return field_vector[position]

    def tokenizer(self, file_path):
        """
        generator of the OF file tokens
        """
        i = 1

        with open(file_path, "r", encoding="utf-8") as file:
            for line in file:
                line = line.strip("\r\n")
                line_tokens = re.finditer(
                    r"[^\s{};()/*<>\[\]]+|{|}|;|\(|\)|//|/\*|\*/|\*|<|>|\[|\]|/|\s+",
                    line,
                )
                # print (line_tokens)

                # if i > 40:
                #    break
                for token in line_tokens:
                    yield (i, token.group())

                i += 1

    def token_classifier(self, file_path):
        """
        generator of the classified OF file tokens

        Args:
            file_path (string): file path
        """
        current_line = 1
        j = 0
        dict_depth = 0
        multi_line_comment = False
        single_line_comment = False
        name = None
        path = "NONE"
        vector_data = False
        list_data = False
        type_data = False
        numerosity = False
        hold_path = ""

        for line, token in self.tokenizer(file_path):
            # print("line: ", line, "current_line: ", current_line, "token: ", token)
            j += 1
            token_dict = {}
            # new line
            if line != current_line:
                # output_file.write(line_string + "\n")
                # line_string = ""
                current_line = line
                # comments
                if multi_line_comment:
                    path = "COMMENT"
                elif single_line_comment:
                    path = hold_path
                    single_line_comment = False

            # identify COMMENTS tokens
            if token == "/*":
                multi_line_comment = True
                if path != "COMMENT":
                    hold_path = path
                path = "COMMENT"
            if token == "*/":
                multi_line_comment = False
            if token == "//" and not multi_line_comment:
                single_line_comment = True
                if path != "COMMENT":
                    hold_path = path
                path = "COMMENT"

            if path == "COMMENT":
                token_type = "COMMENT"
                token_dict["line"] = line
                token_dict["id"] = j
                token_dict["token"] = token
                token_dict["path"] = path
                token_dict["dict_depth"] = dict_depth
                token_dict["type"] = token_type
                if not multi_line_comment and not single_line_comment:
                    path = hold_path
                # print (token_dict)
                # line_string += token
                yield token_dict
                continue

            elif token.isspace():
                token_type = "SPACE"

            elif token == "{":
                token_type = "CURLY_OPEN"
                dict_depth += 1
                if dict_depth == 1:
                    path = name
                    name = None
            elif token == "}":
                token_type = "CURLY_CLOSE"
                dict_depth -= 1
                name = None

            elif token == ";":
                token_type = "FIELD_END"
                name = None
            else:
                if not name:
                    name = token
                    token_type = "FIELD_NAME"
                    if dict_depth == 0:
                        path = name
                    else:
                        path_levels = path.split(".")
                        path = ".".join(path_levels[:dict_depth]) + "." + name
                else:
                    if token == "(":  # list start
                        token_type = "ROUND_OPEN"
                        list_data = True
                    elif token == ")":  # list start
                        token_type = "ROUND_CLOSE"
                        list_data = False
                    elif token == "[":  # list start
                        token_type = "SQUARE_OPEN"
                        vector_data = True
                    elif token == "]":  # list start
                        token_type = "SQUARE_CLOSE"
                        vector_data = False
                    elif token == "<":  # list start
                        token_type = "ANGLE_OPEN"
                        type_data = True
                    elif token == ">":  # list start
                        token_type = "ANGLE_CLOSE"
                        type_data = False
                        numerosity = True
                    elif vector_data:
                        token_type = "VECTOR_DATA"
                    elif list_data:
                        token_type = "LIST_DATA"
                    elif type_data:
                        token_type = "TYPE_DATA"
                    elif numerosity:
                        token_type = "NUMEROSITY"
                        numerosity = False
                    else:
                        token_type = "FIELD_VALUE_GENERIC"

            token_dict["line"] = line
            token_dict["id"] = j
            token_dict["token"] = token
            token_dict["path"] = path
            token_dict["dict_depth"] = dict_depth
            token_dict["type"] = token_type
            yield token_dict

    def convert_phi_to_volumetric(self, out_path):
        """
        This function copy the phi values from the last time folder of the
        of simulation and it copies in another folder with the mass flow rates
        converted into volumetric flow rates.

        Args:
            out_path (string): output file path
        """

        self.write_of_tokens(out_path, self.phi_tokens_to_volumetric())

        return

    def write_of_tokens(self, out_path, token_generator):
        """
        this function writes the tokens in a file

        Args:
            out_path (string): output file path
            tokens (generator): generator of the tokens
        """

        with open(out_path, "w", encoding="utf-8") as output_file:
            line_string = ""
            current_line = 1
            for token_dict in token_generator:
                # new line
                if token_dict["line"] != current_line:
                    output_file.write(line_string + "\n")
                    line_string = ""
                    current_line = token_dict["line"]
                line_string += token_dict["token"]

            output_file.write(line_string + "\n")

        return

    def phi_tokens_to_volumetric(self):
        """
        Generator that converts the phi tokens from mass flow rates to
        volumetric flow rates
        """

        count_vec = 0
        count_internal = 0
        count_boundary = 0
        internal_check = 0
        found_internal = False
        found_boundary = False
        name_boundary = None
        single_val = False
        boundary_check = 0
        i_density = 1 / self.density_kg_m3

        file_path = self.path + "/" + str(self.last_time) + "/phi"

        for t_dict in self.token_classifier(file_path):
            # convert dimensions vector
            if t_dict["path"] == "dimensions" and t_dict["type"] == "VECTOR_DATA":
                count_vec += 1
                if count_vec == 1:
                    t_dict["token"] = str(0)
                if count_vec == 2:
                    t_dict["token"] = str(3)

            if t_dict["path"] == "internalField" and t_dict["type"] == "NUMEROSITY":
                internal_check = int(t_dict["token"])
                found_internal = True

            if t_dict["path"] == "internalField" and t_dict["type"] == "LIST_DATA":
                count_internal += 1
                new_phi = float(t_dict["token"]) * i_density
                t_dict["token"] = f"{new_phi:.11e}"

            if found_internal and t_dict["type"] == "ROUND_CLOSE":
                if internal_check != count_internal:
                    raise ValueError(
                        "ERROR with the conversion of the \
                                     internal phi values"
                    )
                found_internal = False

            if (
                "boundaryField" in t_dict["path"]
                and "value" in t_dict["path"]
                and t_dict["path"] != name_boundary
            ):
                found_boundary = True
                name_boundary = t_dict["path"]

            if found_boundary and t_dict["token"] == "uniform":
                single_val = True

            if found_boundary and t_dict["type"] == "NUMEROSITY":
                boundary_check = int(t_dict["token"])

            if found_boundary and t_dict["type"] == "LIST_DATA":
                count_boundary += 1
                new_phi = float(t_dict["token"]) * i_density
                t_dict["token"] = f"{new_phi:.11e}"

            if found_boundary and t_dict["type"] == "ROUND_CLOSE":
                if boundary_check != count_boundary:
                    raise ValueError(
                        "ERROR with the conversion of the \
                                      boundary phi values"
                    )
                count_boundary = 0

            if found_boundary and single_val and check_float(t_dict["token"]):
                if t_dict["token"] != "0":
                    new_phi = float(t_dict["token"]) * i_density
                    t_dict["token"] = f"{new_phi:.11e}"
                single_val = False

            if found_boundary and t_dict["type"] == "CURLY_CLOSE":
                found_boundary = False

            yield t_dict

    def scale_mesh_results(self, out_path, inlet_activity, decay_const):
        """
        this function takes the vtk um-mesh as computed by fluned, it scales
        it to a new inlet activity, and writes the scaled results to a new vtk file
        """

        print("writing vtk file for circuit calculation ... ")

        self.scaled_vtk_file_path = out_path

        mesh = pv.read(self.vtk_file_path)

        decay_array = mesh.cell_data["Td"]
        rr_array = mesh.cell_data["Ta"]

        del mesh.cell_data["Td"]
        del mesh.cell_data["Ta"]
        del mesh.cell_data["T"]
        # del mesh.cell_data['CellID']
        del mesh.point_data["T"]
        del mesh.point_data["Ta"]
        del mesh.point_data["Td"]

        decay_array = inlet_activity * decay_array / self.inlet_td_conc_atoms_m3[-1]
        rr_array = decay_const * rr_array

        mesh.cell_data["average_vol_activity_bq_m3_decay"] = decay_array
        mesh.cell_data["average_vol_activity_bq_m3_rr"] = rr_array
        mesh.cell_data["average_vol_activity_bq_m3"] = rr_array + decay_array

        mesh.save(out_path)

        return

    def calculate_cartesian_sampling_coordinates(self, sampling_res_cm):
        """
        this function returns a list of dictionaries with the info relative to
        the sampling coordinates
        """

        vtk_boundaries = self.vtk_dimensions

        # convert to cm
        x_bounds = [
            math.floor(vtk_boundaries[0] * 100),
            math.ceil(vtk_boundaries[1] * 100),
        ]
        y_bounds = [
            math.floor(vtk_boundaries[2] * 100),
            math.ceil(vtk_boundaries[3] * 100),
        ]
        z_bounds = [
            math.floor(vtk_boundaries[4] * 100),
            math.ceil(vtk_boundaries[5] * 100),
        ]

        x_ints = math.ceil((x_bounds[1] - x_bounds[0]) / (sampling_res_cm))
        y_ints = math.ceil((y_bounds[1] - y_bounds[0]) / (sampling_res_cm))
        z_ints = math.ceil((z_bounds[1] - z_bounds[0]) / (sampling_res_cm))

        x_nodes = [(x_bounds[0] + sampling_res_cm * i) for i in range(x_ints + 1)]
        y_nodes = [(y_bounds[0] + sampling_res_cm * i) for i in range(y_ints + 1)]
        z_nodes = [(z_bounds[0] + sampling_res_cm * i) for i in range(z_ints + 1)]

        xyz_nodes = [x_nodes, y_nodes, z_nodes]
        sample_coordinates_m = []

        voxel_list = []
        voxel_id = 1

        for i in range(x_ints):
            # xVoxelNodes = [x_nodes[i], x_nodes[i+1]]
            x_voxel_center = (x_nodes[i + 1] + x_nodes[i]) / 2
            for j in range(y_ints):
                # yVoxelNodes = [y_nodes[j], y_nodes[j+1]]
                y_voxel_center = (y_nodes[j + 1] + y_nodes[j]) / 2
                for k in range(z_ints):
                    # zVoxelNodes = [z_nodes[k], z_nodes[k+1]]
                    z_voxel_center = (z_nodes[k + 1] + z_nodes[k]) / 2
                    new_dict = {}
                    new_dict["id"] = voxel_id
                    coords = [x_voxel_center, y_voxel_center, z_voxel_center]
                    new_dict["cent_coords_cm"] = coords
                    coords_m = [coord / 100 for coord in coords]
                    new_dict["cent_coords_m"] = coords_m
                    sample_coordinates_m.append(coords_m)
                    new_dict["volume_cm3"] = sampling_res_cm**3
                    voxel_list.append(new_dict)
                    voxel_id += 1

        self.x_ints = x_ints
        self.y_ints = y_ints
        self.z_ints = z_ints

        self.x_nodes = x_nodes
        self.y_nodes = y_nodes
        self.z_nodes = z_nodes
        self.cartesian_voxel_list = voxel_list
        self.cartesian_sample_coordinates_m = sample_coordinates_m
        self.cartesian_sampling_precision_cm = sampling_res_cm
        self.cartesian_voxel_volume_cm3 = sampling_res_cm**3

        return

    def sample_source_to_cartesian_mesh(self, dataset):
        """
        this function takes the cartesian sampling coordinates computed before and
        use them to:
        1. sample the concentrations from the fluned simulation um mesh
        2. calculate the total sampled values
        3. the sampled emission rates are adjusted so the total emission remains the same

        NB the generated sampled values are in cm scale as this process is done to
        generate radiation sources either for openMC or MCNP
        """

        sampled_cartesian_concentrations = sample_vtk(
            self.vtk_file_path, dataset, self.cartesian_sample_coordinates_m
        )

        voxel_volume = self.cartesian_sampling_precision_cm**3

        for voxel, concentration in zip(
            self.cartesian_voxel_list, sampled_cartesian_concentrations
        ):
            voxel["emission_rate_cm3"] = (
                concentration
                * voxel_volume
                * self.decay_constant
                * self.branching_ratio
                * 1e-06
            )  # atoms per m3 to cm3

        self.raw_sampled_total_emission_rate = sum(
            [voxel["emission_rate_cm3"] for voxel in self.cartesian_voxel_list]
        )

        ratio_vtk_sampling = (
            self.total_isotope_emission_rate / self.raw_sampled_total_emission_rate
        )

        for voxel in self.cartesian_voxel_list:
            voxel["normalized_emission_rate_cm3"] = (
                voxel["emission_rate_cm3"] * ratio_vtk_sampling
            )

        self.normalized_sampled_total_emission_rate = sum(
            [
                voxel["normalized_emission_rate_cm3"]
                for voxel in self.cartesian_voxel_list
            ]
        )

        return

    def generate_system_files(self, time_treatment, fv_scheme):
        """
        this function creates the files needed in the system folder for an
        openFOAM simulation - in later development it will apply the case
        parameters
        """

        # for patch in self.patches.values():
        #     print(patch)
        #     print(patch.face_type)

        control_dict_text_transient = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application    FLUNED-solver;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         100;

deltaT          0.1;

writeControl    timeStep;

writeInterval   100;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

adjustTimeStep  yes;

maxCo           0.5;

includeDecayScalar true;




functions
{

"""

        control_dict_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application    FLUNED-solver;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         1000;

deltaT          1;

writeControl    timeStep;

writeInterval   100;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

includeDecayScalar true;




functions
{

"""
        vol_calc_text = """
    volumeCalc
    {{
        type            writeCellVolumes;
        libs            ("libfieldFunctionObjects.so");
	    select      all;

        writeFields     false;
        writeControl {};

    }}

 """

        vol_flow_text = """
    {}
    {{

        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        patch   {};
        fields  (phi);

        operation sum;
        select  patch;
        name        $patch;

        writeFields     false;
        writeControl {};

    }}

 """
        vol_tx_flow_text = """
    {}
    {{

        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        patch   {};
        fields  ({});
        weightField phi;

        operation sum;
        select  patch;
        name        $patch;

        writeFields     false;
        writeControl {};

    }}

 """

        vol_tx_sum_text = """
    {}
    {{

        type            volFieldValue;
        libs            ("libfieldFunctionObjects.so");
        fields  ({});

        operation       volAverage;
	    select      all;

        writeFields     false;
        writeControl {};


    }}

 """
        system_folder = os.path.join(self.path, "system")
        control_dict_path = os.path.join(system_folder, "controlDict")

        with open(control_dict_path, "w", encoding="utf-8") as fw:
            if time_treatment == "steadystate":
                fw.write(control_dict_text)
                write_control = "outputTime"
            elif time_treatment == "transient":
                fw.write(control_dict_text_transient)
                write_control = "timeStep"

            fw.write(vol_calc_text.format(write_control))

            fw.write(vol_tx_sum_text.format("volTSum", "T", write_control))
            fw.write(vol_tx_sum_text.format("volTaSum", "Ta", write_control))
            fw.write(vol_tx_sum_text.format("volTdSum", "Td", write_control))

            for face in self.patches.values():
                if face.face_type in ["inlet", "outlet"]:
                    fw.write(
                        vol_flow_text.format(
                            "volFlow-" + face.face_id, face.face_id, write_control
                        )
                    )

                    fw.write(
                        vol_tx_flow_text.format(
                            "volTFlow-" + face.face_id,
                            face.face_id,
                            "T",
                            write_control,
                        )
                    )
                    fw.write(
                        vol_tx_flow_text.format(
                            "volTrFlow-" + face.face_id,
                            face.face_id,
                            "Tr",
                            write_control,
                        )
                    )
                    fw.write(
                        vol_tx_flow_text.format(
                            "volTdFlow-" + face.face_id,
                            face.face_id,
                            "Td",
                            write_control,
                        )
                    )
                    fw.write(
                        vol_tx_flow_text.format(
                            "volTaFlow-" + face.face_id,
                            face.face_id,
                            "Ta",
                            write_control,
                        )
                    )

            fw.write("}")

        fv_scheme_stable = """
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}}

ddtSchemes
{{
    default         {};
}}

gradSchemes
{{
    default         cellLimited Gauss linear 1;
}}

divSchemes
{{
    default         none;
    div(phi,T)      Gauss  upwind;
    div(phi,Ta)      Gauss upwind;
    div(phi,Td)      Gauss upwind;
    div(phi,Tr)      Gauss upwind;
}}

laplacianSchemes
{{
    default         none;
    laplacian(DT,T) Gauss linear corrected;
    laplacian(Dturbulent,T) Gauss linear corrected;
    laplacian(DT,Tr) Gauss linear corrected;
    laplacian(Dturbulent,Tr) Gauss linear corrected;
    laplacian(DT,Ta) Gauss linear corrected;
    laplacian(Dturbulent,Ta) Gauss linear corrected;
    laplacian(DT,Td) Gauss linear corrected;
    laplacian(Dturbulent,Td) Gauss linear corrected;
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         limited 0.5;
}}

"""

        fv_scheme_accurate = """
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}}

ddtSchemes
{{
    default         {};
}}

gradSchemes
{{
    default         cellLimited Gauss linear 0.5;
}}

divSchemes
{{
    default         none;
    div(phi,T)      Gauss linearUpwind default;
    div(phi,Ta)      Gauss linearUpwind default;
    div(phi,Td)      Gauss linearUpwind default;
    div(phi,Tr)      Gauss linearUpwind default;
}}

laplacianSchemes
{{
    default         none;
    laplacian(DT,T) Gauss linear corrected;
    laplacian(Dturbulent,T) Gauss linear corrected;
    laplacian(DT,Tr) Gauss linear corrected;
    laplacian(Dturbulent,Tr) Gauss linear corrected;
    laplacian(DT,Ta) Gauss linear corrected;
    laplacian(Dturbulent,Ta) Gauss linear corrected;
    laplacian(DT,Td) Gauss linear corrected;
    laplacian(Dturbulent,Td) Gauss linear corrected;
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         limited 1;
}}

"""

        schemes_path = os.path.join(system_folder, "fvSchemes")

        with open(schemes_path, "w", encoding="utf-8") as fw:
            # select the fv scheme
            if fv_scheme == "accurate":
                fv_scheme_text = fv_scheme_accurate
            elif fv_scheme == "stable":
                fv_scheme_text = fv_scheme_stable

            # select the time treatment
            elif time_treatment == "transient":
                fw.write(fv_scheme_text.format("Euler"))
            if time_treatment == "steadystate":
                fw.write(fv_scheme_text.format("steadyState"))
            elif time_treatment == "transient":
                fw.write(fv_scheme_text.format("Euler"))

        fv_solution_text = """
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


solvers
{{
    T
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-07;
        relTol          {};
    }}
    Tr
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-07;
        relTol          {};
    }}
    Ta
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-07;
        relTol          {};
    }}
    Td
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-07;
        relTol          {};
    }}
}}

SIMPLE
{{
    nNonOrthogonalCorrectors 0;
    consistent      yes;

    residualControl
    {{

        T               1e-6;
        Tr              1e-6;
        Ta              1e-6;
        Td              1e-6;

    }}


}}

"""
        solution_path = os.path.join(system_folder, "fvSolution")

        with open(solution_path, "w", encoding="utf-8") as fw:
            if time_treatment == "steadystate":
                fw.write(fv_solution_text.format(0.01, 0.01, 0.01, 0.01))
            elif time_treatment == "transient":
                fw.write(fv_solution_text.format(0, 0, 0, 0))

        parallel_dict_text = """
FoamFile
{
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 4;

method          scotch;

// *********************************************************************** //
"""
        parallel_dict_path = os.path.join(system_folder, "decomposeParDict")
        with open(parallel_dict_path, "w", encoding="utf-8") as fw:
            fw.write(parallel_dict_text)

        return

    def generate_constant_file(
        self, isotope, molecular_diffusion, decay_constant, schmidt_number
    ):
        """
        This function creates the files in the constant folder
        """

        transport_prop_text = """
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      transportProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


transportModel  Newtonian;

isotope          {};

DT             DT [0 2 -1 0 0 0 0] {};

lambda         lambda [0 0 -1 0 0 0 0] {};

Sct            Sct [ 0 0 0 0 0 0 0 ] {};


// ************************************************************************ //
"""

        constant_folder = os.path.join(self.path, "constant")
        transport_prop_path = os.path.join(constant_folder, "transportProperties")

        with open(transport_prop_path, "w", encoding="utf-8") as fw:
            fw.write(
                transport_prop_text.format(
                    isotope,
                    molecular_diffusion,
                    decay_constant,
                    schmidt_number,
                )
            )

        return

    def generate_zero_t(self, inlet_conc):
        """
        this function generate the T file at t=0
        """

        zero_folder = os.path.join(self.path, "0")
        zero_t_path = os.path.join(zero_folder, "T")

        t_header_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      T;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}

"""

        with open(zero_t_path, "w", encoding="utf-8") as fw:
            fw.write(t_header_text)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                if face.face_type == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "inlet":
                    fw.write("        type            fixedValue;\n")
                    vals_string = "        value         uniform {};\n"
                    fw.write(vals_string.format(inlet_conc))
                else:
                    raise ValueError(
                        "ERROR: patch type not recognized - ", face.face_type
                    )

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_zero_ta(self):
        """
        this function generate the Ta file at t=0
        """

        zero_folder = os.path.join(self.path, "0")
        zero_ta_path = os.path.join(zero_folder, "Ta")

        ta_header_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Ta;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}
"""

        with open(zero_ta_path, "w", encoding="utf-8") as fw:
            fw.write(ta_header_text)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                if face.face_type == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "inlet":
                    fw.write("        type            fixedValue;\n")
                    vals_string = "        value          uniform 0;\n"
                    fw.write(vals_string)
                else:
                    raise ValueError(
                        "ERROR: patch type not recognized - ", face.face_type
                    )

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_zero_tr(self, time_treatment):
        """
        this function generate the Tr file at t=0
        """

        zero_folder = os.path.join(self.path, "0")
        zero_tr_path = os.path.join(zero_folder, "Tr")

        tr_header_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Tr;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}
"""

        with open(zero_tr_path, "w", encoding="utf-8") as fw:
            fw.write(tr_header_text)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                if face.face_type == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "inlet":
                    fw.write("        type            fixedValue;\n")
                    if time_treatment == "steadystate":
                        vals_string = "        value          uniform 0;\n"
                    elif time_treatment == "transient":
                        vals_string = "        value          uniform 1;\n"
                    fw.write(vals_string)
                else:
                    raise ValueError(
                        "ERROR: patch type not recognized - ", face.face_type
                    )

                fw.write("    }\n")

            fw.write(end_text)

        return

    def generate_zero_td(self, inlet_conc):
        """
        this function generate the Td file at t=0
        """

        zero_folder = os.path.join(self.path, "0")
        zero_td_path = os.path.join(zero_folder, "Td")

        td_header_text = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Td;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

"""
        end_text = """

}

"""

        with open(zero_td_path, "w", encoding="utf-8") as fw:
            fw.write(td_header_text)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                if face.face_type == "wall":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "outlet":
                    fw.write("        type            zeroGradient;\n")
                elif face.face_type == "inlet":
                    fw.write("        type            fixedValue;\n")
                    vals_string = "        value         uniform {};\n"
                    fw.write(vals_string.format(inlet_conc))
                else:
                    raise ValueError(
                        "ERROR: patch type not recognized - ", face.face_type
                    )

                fw.write("    }\n")

            fw.write(end_text)

        return

    def read_volumes(self):
        """
        this function reads the volumes from the V file located in the
        zero folder
        """

        # common patterns
        internalBlockPat = re.compile(
            r"internalField.*?\((.{1,}?)\)", re.MULTILINE | re.DOTALL
        )

        v_file_path = os.path.join(self.path, "0")
        v_files = [f for f in os.listdir(v_file_path) if re.match(r"V(c)?(\..*)?$", f)]
        if not v_files:
            raise FileNotFoundError("No V or Vc files found")

        v_file = os.path.join(v_file_path, v_files[0])

        text = open_utf8_or_gzip(v_file)
        numInternalBlocks = internalBlockPat.findall(text)
        internalVolumes = numInternalBlocks[0].split("\n")[1:-1]
        self.volumes = np.array([float(val) for val in internalVolumes])

        return

    def read_centroids(self):
        """this function reads the centroids from the 0 folder"""

        # common patterns
        internalBlockPat = re.compile(
            "internalField.*?\((.{1,}?)\n\\s*\)", re.MULTILINE | re.DOTALL
        )

        cFile = os.path.join(self.path, "0", "C")
        try:
            inpFile = open(cFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            raise FileNotFoundError("couldn't open  C file")
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalCentroids = numInternalBlocks[0].split("\n")[1:]
            internalCentroids = [val.strip("()") for val in internalCentroids]
            self.centroids = np.array(
                [[float(val) for val in v.split()] for v in internalCentroids]
            )

        return

    def read_velocities(self):
        """this function reads the velocity from the U file located in the
        zero folder"""

        print("reading velocity values...")

        # common patterns
        internalBlockPat = re.compile(
            r"internalField.*?\((.{1,}?)\n\\s*\)", re.MULTILINE | re.DOTALL
        )

        uFile = os.path.join(self.path, "0", "U")
        try:
            inpFile = open(uFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            raise FileNotFoundError("couldn't open  U file")
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalVels = numInternalBlocks[0].split("\n")[1:]
            internalVels = [val.strip("()") for val in internalVels]

            self.velocities = np.array(
                [np.array([float(val) for val in v.split()]) for v in internalVels]
            )

        return

    def read_grad_t(self):
        """this function reads the T gradient"""

        print("reading gradient values...")

        # common patterns
        internalBlockPat = re.compile(
            r"internalField.*?\((.{1,}?)\n\\s*\)", re.MULTILINE | re.DOTALL
        )

        gradFile = os.path.join(self.last_time, "grad(T)")
        try:
            inpFile = open(gradFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            raise FileNotFoundError("couldn't open grad file")
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalVels = numInternalBlocks[0].split("\n")[1:]
            internalVels = [val.strip("()") for val in internalVels]
            self.gradients = np.array(
                [np.array([float(val) for val in v.split()]) for v in internalVels]
            )

        return

    def assign_activation_rates(
        self,
        activation_file,
        activation_const,
        activation_dataset,
        activation_dataset_error,
        activation_normalization,
    ):
        """this file create the source file in the zero folder.
        There are three modes:
        1) no activation, then the file contains only zeros.
        2) constant activation (with input value)
        3) with a source file
        """

        if activation_file == "":
            if activation_const == 0:
                # case with no rad source
                activ_sources = [0 for _ in range(self.n_internal_cells)]
            else:
                # case with constant value source
                self.read_volumes()
                activ_sources = [activation_const * vol for vol in self.volumes]

        else:
            # case with source file
            launch_centroid_func_object(self.path)
            self.read_volumes()
            self.read_centroids()

            # 1.sample the activation file
            print("Sampling Reaction Rate file ... ")
            sampledRates = sample_vtk(
                activation_file, activation_dataset, self.centroids
            )
            # remove possible small negative values
            sampledRates = [val if val > 0 else 0 for val in sampledRates]

            # 1.1 if present sample the vtk file to get the error array
            if activation_dataset_error != "":
                print("Sampling Reaction Rate MCNP errors  ... ")
                self.sampled_rr_error = sample_vtk(
                    activation_file, activation_dataset_error, self.centroids
                )

            # 2.use the activation const as a factor
            if activation_const == 0:
                factor = 1
            else:
                factor = activation_const
            activ_sources = [factor * rate for rate in sampledRates]

            # 3.apply a normalization factor if provided
            if activation_normalization != 0:
                vec = [rate * vol for rate, vol in zip(activ_sources, self.volumes)]

                totalSampled = sum(vec)
                print("total sampled atoms/s")
                print(totalSampled)
                normFactor = activation_normalization / totalSampled

                activ_sources = [rate * normFactor for rate in activ_sources]

                nVec = [rate * vol for rate, vol in zip(activ_sources, self.volumes)]

                print("new total sampled atoms/s")
                print(sum(nVec))

        self.reaction_rates = activ_sources

        return None

    def generate_source_file(self, activation_dataset_error=""):
        """
        generate the source file using the dataset reaction rate data
        """

        zeroFolder = os.path.join(self.path, "0")
        zeroSourcePath = os.path.join(zeroFolder, "Source")

        sHeaderText = """
        FoamFile
        {
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0";
            object      Source;
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        dimensions      [0 0 -1 0 0 0 0];

        internalField   nonuniform List<scalar>
        """
        boundaryText = """
        boundaryField
        {

        """
        closerText = """

        }


        // *********************************************************************** //
        """

        with open(zeroSourcePath, "w") as fw:
            fw.write(sHeaderText)
            fw.write("{:d}\n".format(self.n_internal_cells))
            fw.write("(\n")
            for val in self.reaction_rates:
                fw.write("{:e}\n".format(val))
            fw.write(")\n;\n\n")

            fw.write(boundaryText)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                fw.write("        type            fixedValue;\n")
                valString = "        value          uniform 0;\n"
                fw.write(valString)
                fw.write("    }\n")

            fw.write(closerText)

        if activation_dataset_error != "":
            zeroSourceErrorPath = os.path.join(zeroFolder, "SourceError")

            eHeaderText = """
        FoamFile
        {
            version     2.0;
            format      ascii;
            class       volScalarField;
            location    "0";
            object      SourceError;
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        dimensions      [0 0 -1 0 0 0 0];

        internalField   nonuniform List<scalar>
        """
            with open(zeroSourceErrorPath, "w") as fw:
                fw.write(eHeaderText)
                fw.write("{:d}\n".format(self.n_internal_cells))
                fw.write("(\n")
                for val in self.sampled_rr_error:
                    fw.write("{:e}\n".format(val))
                fw.write(")\n;\n\n")

                fw.write(boundaryText)

                for face in self.patches.values():
                    fw.write("    " + face.face_id + "\n    {\n")
                    fw.write("        type            fixedValue;\n")
                    valString = "        value          uniform 0;\n"
                    fw.write(valString)
                    fw.write("    }\n")

                fw.write(closerText)

        return

    def generate_tr_source_file(self, time_treatment):
        """
        this file create the source file for the time residency fictiotious
        scalar in the zero folder.
        There are two modes:
        1) steady state: the source is 1 in every cell
        2) transient the source is 0 in every cell
        """

        zeroFolder = os.path.join(self.path, "0")
        zeroSourcePath = os.path.join(zeroFolder, "TrSource")

        sHeaderText = """
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      TrSource;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

"""
        intFieldText = """

internalField   uniform {};

        """
        boundaryText = """
boundaryField
{

"""
        closerText = """

}


// *********************************************************************** //
"""

        with open(zeroSourcePath, "w") as fw:
            fw.write(sHeaderText)

            if time_treatment == "steadystate":
                fw.write(intFieldText.format(1))
            elif time_treatment == "transient":
                fw.write(intFieldText.format(0))

            fw.write(boundaryText)

            for face in self.patches.values():
                fw.write("    " + face.face_id + "\n    {\n")
                fw.write("        type            fixedValue;\n")
                valString = "        value          uniform 0;\n"
                fw.write(valString)
                fw.write("    }\n")

            fw.write(closerText)

        return

    def assign_isotope(self):
        """
        using the decay constant this function understand if we are
        considering N-16, N-17 or O-19 and assign the spectrum
        accordingly. If it is not possible a dummy spectrum is
        assigned

        this will be changed by embedding the data in the simulation
        """

        N16_decay_constant = 0.09721559
        N17_decay_constant = 0.1661825
        O19_decay_constant = 0.02578672546
        F20_decay_constant = 0.062093271

        if self.isotope == "custom":
            if math.isclose(self.decay_constant, N16_decay_constant, rel_tol=1e-3):
                self.isotope = "n16"
            elif math.isclose(self.decay_constant, N17_decay_constant, rel_tol=1e-3):
                self.isotope = "n17"
            elif math.isclose(self.decay_constant, O19_decay_constant, rel_tol=1e-3):
                self.isotope = "o19"
            elif math.isclose(self.decay_constant, F20_decay_constant, rel_tol=1e-3):
                self.isotope = "f20"
            else:
                raise ValueError("ERROR, isotope not recognized")

        isotope_database = load_isotopes()
        if self.isotope.lower() not in isotope_database:
            raise ValueError("ERROR isotope not found in the database")
        isotope_data = isotope_database[self.isotope.lower()]
        self.e_bins = isotope_data.e_bins
        self.p_bins = isotope_data.p_bins
        self.e_lines = isotope_data.e_lines
        self.p_lines = isotope_data.p_lines
        self.branching_ratio = isotope_data.branching_ratio
        self.particle_type = isotope_data.emitting_particle

        return

    def parse_constants_file(self):
        """this function parses the constant properties to get the decay
        variable and the others"""

        # common patterns
        dtPat = re.compile("DT\s*DT.*")
        isotope_pat = re.compile("isotope.*")
        lambdaPat = re.compile("lambda\s*lambda.*")
        schPat = re.compile(r"Sct\s*Sct.*")

        cFile = os.path.join(self.path, "constant", "transportProperties")
        try:
            inpFile = open(cFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            raise FileNotFoundError("couldn't open transportProperties file")
        with inpFile:
            text = inpFile.read()
            dtLines = dtPat.findall(text)
            lambdaLines = lambdaPat.findall(text)
            schLines = schPat.findall(text)
            isotopeLines = isotope_pat.findall(text)
            if len(isotopeLines) != 0:
                vals = isotopeLines[0].strip(" ;").split()
                val = vals[-1]
                self.isotope = val
            else:
                self.isotope = "custom"

            if len(dtLines) != 0:
                vals = dtLines[0].strip(" ;").split()
                val = vals[-1]
                self.molecular_diffusion = float(val)
            else:
                self.molecular_diffusion = 0

            if len(lambdaLines) != 0:
                vals = lambdaLines[0].strip(" ;").split()
                val = vals[-1]
                self.decay_constant = float(val)
            else:
                self.decay_constant = 0

            if len(schLines) != 0:
                vals = schLines[0].strip(" ;").split()
                val = vals[-1]
                self.schmidt_number = float(val)
            else:
                self.schmidt_number = 0

        return

    def get_time_treatment(self):
        """
        this function parse the fvSchemes file to check if the simulation
        is steady state or transient
        """

        # common patterns
        ddtPat = re.compile("ddtSchemes.*?\{.{1,}?\}", re.MULTILINE | re.DOTALL)

        cFile = os.path.join(self.path, "system", "fvSchemes")
        try:
            inpFile = open(cFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            raise FileNotFoundError("couldn't open fvSchemes file")
        with inpFile:
            text = inpFile.read()
            ddtText = ddtPat.findall(text)[0]
            if "euler" in ddtText.lower():
                self.time_treatment = "transient"
            elif "steadystate" in ddtText.lower():
                self.time_treatment = "steadystate"
            else:
                raise ValueError("couldn't open recognise the time scheme")

        return

    def read_t(self):
        """
        this function reads the T scalar
        """

        print("reading scalar values...")

        # find the last folder
        # folderItms = os.listdir(self.fluned_path)
        # folderTimes=[int(itm) for itm in folderItms if checkInt(itm) == True]
        # lastTime = max(folderTimes)

        # common patterns
        internalBlockPat = re.compile(
            r"internalField.*?\((.{1,}?)\)", re.MULTILINE | re.DOTALL
        )

        tFile = os.path.join(self.last_time, "T")

        try:
            inpFile = open(tFile, "r", encoding="utf8", errors="ignore")
        except IOError:
            raise FileNotFoundError("couldn't open T file")
        with inpFile:
            text = inpFile.read()
            numInternalBlocks = internalBlockPat.findall(text)
            internalScalar = numInternalBlocks[0].split("\n")[1:-1]

            self.t_scalar = np.zeros(self.n_internal_cells)
            for i in range(self.n_internal_cells):
                self.t_scalar[i] = float(internalScalar[i])

        return

    def create_results_folder(self):
        """
        this function creates the results folder
        """

        resFolder = os.path.join(self.path, "RESULTS")

        dir_check = os.path.isdir(resFolder)
        if not dir_check:
            os.mkdir(resFolder)

        self.results_folder = resFolder

        return

    def write_sampled_cartesian_source_vtk(self):
        """
        this function write the sampled source in a vtk file
        """

        out_vtk_file_path = os.path.join(
            self.results_folder, "cartesian_sampled_source.vtk"
        )

        dims = (self.x_ints + 1, self.y_ints + 1, self.z_ints + 1)

        spacing = (
            (self.x_nodes[-1] - self.x_nodes[0]) / self.x_ints,
            (self.y_nodes[-1] - self.y_nodes[0]) / self.y_ints,
            (self.z_nodes[-1] - self.z_nodes[0]) / self.z_ints,
        )

        origin = (self.x_nodes[0], self.y_nodes[0], self.z_nodes[0])

        raw = np.asarray(
            [vox["normalized_emission_rate_cm3"] for vox in self.cartesian_voxel_list],
            dtype=float,
        )
        reordered = raw.reshape(
            (self.x_ints, self.y_ints, self.z_ints), order="C"
        ).ravel(order="F")

        dataset_name_output = "emission_rate_cm3"

        write_cartesian_vtk(
            out_vtk_file_path, dims, spacing, origin, reordered, dataset_name_output
        )

        return

    def write_cdgs(self):
        """
        this function write the sample CDGS file
        """

        cdgsFile = os.path.join(self.results_folder, "cartesian_sampled_source.cdgs")

        with open(cdgsFile, "w") as fw:
            fw.write("num_meshes 1\n")
            fw.write(
                "global_source {:e}\n".format(
                    self.normalized_sampled_total_emission_rate
                )
            )
            fw.write("mesh_id 1\n")
            fw.write("cdgs from vtk\n")
            fw.write("Cooling_time 0.0\n")
            fw.write(
                "total_source {:e}\n".format(
                    self.normalized_sampled_total_emission_rate
                )
            )
            fw.write("energy_type {}\n".format("bins"))
            fw.write("energy_boundaries {:d}\n".format(len(self.e_bins)))
            # WRITE SPECTRUM BINS
            specString = formatValues(self.e_bins)
            fw.write(specString)

            fw.write("mesh_type rec\n")
            fw.write(
                "mesh_boundaries {:d} {:d} {:d}\n".format(
                    self.x_ints + 1, self.y_ints + 1, self.z_ints + 1
                )
            )
            fw.write("0.000000e+00  0.000000e+00  0.000000e+00\n")
            fw.write("1.000000e+00  0.000000e+00  0.000000e+00\n")
            fw.write("0.000000e+00  1.000000e+00  0.000000e+00\n")
            xString = formatValues(self.x_nodes)
            fw.write(xString)
            yString = formatValues(self.y_nodes)
            fw.write(yString)
            zString = formatValues(self.z_nodes)
            fw.write(zString)
            fw.write("source_data\n")

            voxelString1 = "{:d} {:.5e} {:.5e} 1\n"
            voxelString2 = "0 1.0 {:.5e}\n"
            specErrorString = formatValues([0] * (len(self.p_bins)))

            for vox in self.cartesian_voxel_list:
                if vox["normalized_emission_rate_cm3"] > 0:
                    fw.write(
                        voxelString1.format(
                            vox["id"],
                            vox["normalized_emission_rate_cm3"],
                            self.cartesian_voxel_volume_cm3,
                        )
                    )
                    fw.write(voxelString2.format(vox["normalized_emission_rate_cm3"]))
                    emittingSpectrum = [
                        val * vox["normalized_emission_rate_cm3"] for val in self.p_bins
                    ]
                    spectrumString = formatValues(emittingSpectrum)
                    fw.write(spectrumString)
                    fw.write(specErrorString)

            fw.write("end_source_data\n")

        return
