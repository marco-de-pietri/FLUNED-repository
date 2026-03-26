"""
oFoamBase class which is base class for the class manipulating the openfoam simulations
"""

import re
from pathlib import Path

import numpy as np
from foamlib import FoamCase, FoamFieldFile

from .patch_class import PatchClass


def check_float(s):
    """
    check if a string can be converted to float
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


class oFoamBase:
    """
    A class to represent a simulation using OpenFOAM.
    some functions are present to parse and generate file

    Attributes:
    -----------
    path : str
        The complete path to the simulation folder.
    """

    def __init__(self, path: str | Path):
        """
        Constructs all the necessary attributes for the oFoamBase object.

        Parameters:
        -----------
        path : str
            The complete path to the simulation folder.
        """
        self.path = Path(path)
        self.case = self.path.name
        self.foamlib_object = FoamCase(self.path)

    def create_case_folder(self):
        """
        Create an OpenFOAM case folder structure.
        """

        case_path = Path(self.path)

        # Create required directory structure
        for subpath in [
            case_path,
            case_path / "0",
            case_path / "constant" / "polyMesh",
            case_path / "system",
        ]:
            subpath.mkdir(parents=True, exist_ok=True)

        # Create case.foam file
        (case_path / "case.foam").touch(exist_ok=True)

    @classmethod
    def is_valid_openfoam_case_directory(cls, folder_path: str | Path) -> bool:
        """
        Check if a path is a completed OpenFOAM simulation (case) directory.
        Returns True if all of the following subdirectories exist:
        - "0" (initial time directory)
        - at least one other numeric directory (time > 0)
        - "constant" (mesh & physical properties)
        - "system" (control dictionaries)
        """
        folder_path = Path(folder_path)
        if not folder_path.is_dir():
            return False

        entries = {item.name for item in folder_path.iterdir() if item.is_dir()}

        has_initial_time = "0" in entries
        has_later_time = any(name.isdigit() and int(name) > 0 for name in entries)
        has_constant = "constant" in entries
        has_system = "system" in entries

        return has_initial_time and has_later_time and has_constant and has_system

    def post_process_openfoam_simulation(self):
        """
        this function is called when we process a finished openfoam simulation
        """

        self.post_process_path = self.path / "postProcessing"
        self.last_time = self.foamlib_object[-1].name

        self.volumetric_flag = self.get_volumetric_flag()

        if not self.volumetric_flag:
            self.density_kg_m3 = self.get_density()
        else:
            self.density_kg_m3 = 1000

        self.patches = self.get_patches()
        self.n_internal_cells = self.get_number_internal_cells()

    def get_number_internal_cells(self) -> int:
        """
        this function create an attribute to the class that specifies the
        number of internal cells. It does so by reading the U file
        """

        u_0 = self.path / "0" / "U"
        u_last = self.path / self.last_time / "U"

        num_cells_0 = 0
        num_cells_last = 0

        if u_0.is_file():
            internal_field = FoamFieldFile(u_0).internal_field
            if internal_field.shape == (3,):
                num_cells_0 = 0
            else:
                num_cells_0 = internal_field.shape[0]

        if u_last.is_file():
            internal_field = FoamFieldFile(u_last).internal_field
            if internal_field.shape == (3,):
                num_cells_last = 0
            else:
                num_cells_last = internal_field.shape[0]

        if num_cells_0 == 0 and num_cells_last == 0:
            raise ValueError("Unable to determine number of internal cells")

        if (num_cells_0 != 0 and num_cells_last != 0) and (
            num_cells_0 != num_cells_last
        ):
            raise ValueError(
                "Error number of internal cells in 0 and last time are different"
            )

        return int(num_cells_0) if num_cells_0 != 0 else int(num_cells_last)

    def get_patches(self):
        """
        this function creates a dictionary with the patches objects
        """

        patches = {}

        patches_list = self.parse_boundary_phi_files()

        for patch in patches_list:
            patches[patch["face_id"]] = PatchClass(patch, self.path)

        return patches

    def parse_boundary_phi_files(self):
        """
        this function examines phi files and returns a list of dictionary with
        the face_name and type
        at the moment it is not optimized for speed
        """

        # print ("parsing faces...")

        faces = []

        poly_mesh_folder = self.path / "constant" / "polyMesh"
        boundary_file = poly_mesh_folder / "boundary"

        try:
            inp_file = open(boundary_file, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open boundary file")
            raise FileNotFoundError(f"Unable to open phi file at '{boundary_file}'")
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

        phi_file_path = self.path / str(self.last_time) / "phi"
        try:
            inp_file = open(phi_file_path, "r", encoding="utf8", errors="ignore")
        except IOError:
            print("couldn't open phi file")
            raise FileNotFoundError(f"Unable to open phi file at '{phi_file_path}'")
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

    def get_volumetric_flag(self):
        """
        this function check the phi values to check if the phis are in the units
        of volumetric or mass flux (m3/s or kg/s)
        """

        phi_dimensions = tuple(self.foamlib_object[self.last_time]["phi"].dimensions)

        if phi_dimensions == (1, 0, -1, 0, 0, 0, 0):
            volumetric_flag = False
        elif phi_dimensions == (0, 3, -1, 0, 0, 0, 0):
            volumetric_flag = True
        else:
            raise ValueError(
                "The dimensions of phi are not recognized: ", phi_dimensions
            )

        return volumetric_flag

    def get_density(self):
        """
        this function gets the density value from the
        constant/transportProperties file
        """

        file_path = self.path / "constant" / "transportProperties"
        dict_path = "rhoRef"

        value = self.query_of_single_value(file_path, dict_path)

        return float(value)

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

        file_path = self.path / str(self.last_time) / "phi"

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

    def read_volumes(self):
        """
        this function reads the volumes from the V file located in the
        zero folder
        """

        time_files = self.foamlib_object["0"]
        if "V" in time_files:
            vol_name = "V"
        elif "Vc" in time_files:
            vol_name = "Vc"
        else:
            raise FileNotFoundError("No V or Vc files found")

        self.volumes = self.foamlib_object["0"][vol_name].internal_field

        return

    def read_centroids(self):
        """Read centroids from the OpenFOAM '0/C' file and store them in self.centroids (N x 3)."""

        self.centroids = self.foamlib_object["0"].cell_centers().internal_field

    def read_velocities(self):
        """this function reads the velocity from the U file located in the
        zero folder"""

        print("reading velocity values...")

        self.velocities = self.foamlib_object["0"]["U"].internal_field

        return

    def get_time_treatment(self):
        """
        Parse the fvSchemes file to determine whether the simulation
        is steady-state or transient.
        """

        default_ddt = self.foamlib_object.fv_schemes["ddtSchemes", "default"].lower()

        if "euler" in default_ddt:
            self.time_treatment = "transient"
        elif "steadystate" in default_ddt:
            self.time_treatment = "steadystate"
        else:
            raise ValueError("could not recognize the time discretization scheme")

        return self.time_treatment
