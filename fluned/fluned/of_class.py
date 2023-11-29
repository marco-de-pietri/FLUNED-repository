"""
class for the OF simulations, this can be used to parse and generate files
"""
import os
import re
import sys
import copy
import numpy as np

def is_float(s):
    """
    check if a string can be converted to float
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

def merge_continue_runs(time_lists,data_lists):
    """
    this function takes time series segments and join them in a single one -
    the overlapping sections are removed
    """

    if len(time_lists) == 0 or len(data_lists) == 0 :
        raise ValueError('Error with the number of post process files')

    if len(time_lists) != len(data_lists) :
        print ("ERROR mismatch in length of data series")
        print ("number of time lists: ", len(time_lists))
        print ("number of data lists: ", len(data_lists))
        raise ValueError('Error with the number of post process files')

    tot_len_time_lists = sum([len(x) for x in time_lists])
    tot_len_data_lists = sum([len(x) for x in data_lists])

    if tot_len_time_lists != tot_len_data_lists:
        print ("ERROR mismatch in length of data series")
        print ("number of time points: ", tot_len_time_lists)
        print ("number of data points: ", tot_len_data_lists)
        raise ValueError('Error with the number of post process data points')


    temp_list = list(zip(time_lists,data_lists))

    #time_lists_sorted = sorted(time_lists,
    #                           key=lambda x:x[0])

    sorted_temp_list = sorted(temp_list, key= lambda x:x[0][0])

    time_lists_sorted,data_lists_sorted = zip(*sorted_temp_list)
    #data_lists_sorted = sorted(data_lists,
    #                           key=lambda x:x[0])

    time_series = time_lists_sorted[0]
    data_series = data_lists_sorted[0]


    if len(time_lists) > 1:
        for i,l in enumerate(time_lists_sorted):
            time_series_temp = copy.deepcopy(time_series)
            if i == 0:
                continue
            appended = False
            for j, val in enumerate(time_series_temp):
                if val >= time_lists_sorted[i][0]:
                    time_series = time_series[0:j]
                    time_series.extend(l)
                    data_series = data_series[0:j]
                    data_series.extend(data_lists_sorted[i])
                    appended = True
                    break

            if not appended:
                time_series.extend(l)
                data_series.extend(data_lists_sorted[i])


    return time_series, data_series

def post_file_array(path_list,name):
    """
    this function extracts a generic array contained in the post
    processing file
    """

    array_list = []

    for f_path in path_list:

        array = []
        index = -1

        try:
            post_file = open(f_path,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open postprocess file")
        with post_file:
            lines = post_file.readlines()
            for line in lines:

                line = line.replace('#','')

                wrds = line.split()

                if len(wrds) == 0:
                    continue

                if 'time' in wrds[0].lower():
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
    folder_itms = sorted(folder_itms, reverse=False,key=lambda x:float(x))

    for fld in folder_itms:

        path_1 = os.path.join(path,fld)
        if os.path.isdir(path_1):
            file_itms = os.listdir(path_1)
            complete_path = os.path.join(path_1,file_itms[0])
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
        self.post_process_path = os.path.join(self.simulation_path,'postProcessing')
        self.face_id = patch_dict['face_id']
        self.face_elements_n = patch_dict['face_elements_n']
        self.face_first_element = patch_dict['face_first_element']
        self.face_type = patch_dict['type']
        self.face_sum_phis = patch_dict['sum_phis']
        self.post_process_flow = []
        self.post_process_time = []
        self.post_process_t_flow = []
        self.post_process_ta_flow = []
        self.post_process_td_flow = []
        self.t_conc_atoms_m3 = []
        self.ta_conc_atoms_m3 = []
        self.td_conc_atoms_m3 = []

    def post_process_face(self):
        """
        this function reads all the post process files and calculates the
        different scalar flows across the face
        """
        self.post_process_flow, self.post_process_time = self.get_post_process_list('volFlow-','sum(phi)')
        self.post_process_t_flow, _ = self.get_post_process_list('volTFlow-','sum(T)')
        self.post_process_ta_flow, _ = self.get_post_process_list('volTaFlow-','sum(Ta)')
        self.post_process_td_flow, _ = self.get_post_process_list('volTdFlow-','sum(Td)')
        self.t_conc_atoms_m3 = [x/y if y!= 0 else 0 for x,y in zip(
                                    self.post_process_t_flow,
                                    self.post_process_flow)]
        self.ta_conc_atoms_m3 = [x/y if y!= 0 else 0 for x,y in zip(
                                    self.post_process_ta_flow,
                                    self.post_process_flow)]
        self.td_conc_atoms_m3 = [x/y if y!= 0 else 0 for x,y in zip(
                                    self.post_process_td_flow,
                                    self.post_process_flow)]


        # check the sign of the flow
        if self.post_process_flow[-1] > 0 and self.face_type != 'outlet':
            raise ValueError('Error with the flow sign')
        elif self.post_process_flow[-1] < 0 and self.face_type != 'inlet':
            raise ValueError('Error with the flow sign')
        elif self.post_process_flow[-1] == 0 and self.face_type != 'wall':
            raise ValueError('Error with the flow sign')
        if self.t_conc_atoms_m3[-1] < 0:
            raise ValueError('Error with the concentration sign')

        print (self.post_process_t_flow)
        print (self.t_conc_atoms_m3)

        return None


    def get_post_process_list(self, dir_prefix, col_name ):
        """
        this function gets the post process flow from the simulation folder
        """

        if self.face_type == 'wall':
            return [0],[0]

        dir_name = dir_prefix + self.face_id
        flow_folders = [itm for itm in os.listdir(self.post_process_path)
                        if dir_name in itm]

        if len(flow_folders) != 1:
            raise ValueError('Error with the number of flow folders')

        #print (flow_folders)

        flow_files = get_post_files(os.path.join(self.post_process_path,
                                                flow_folders[0]))

        #print (flow_files)
        #postDic['areaFile'] = postFileArea(postDic['flowFiles'][0])

        time_lists = post_file_array(flow_files,'Time')
        flow_lists = post_file_array(flow_files,col_name)


        time_list_sorted,flow_list_sorted = merge_continue_runs (
                                        time_lists,
                                        flow_lists,
                                        )



        return flow_list_sorted, time_list_sorted


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
        self.last_time = self.get_last_time()

        self.volumetric_flag = self.get_volumetric_flag()


        if not self.volumetric_flag:
            self.density_kg_m3 = self.get_density()
        else:
            self.density_kg_m3 = 1000

        self.patches = self.get_patches()
        self.reduction_rate = 0


    def post_process_simulation(self):
        """
        this function is called when we process a finished FLUNED simulation
        """
        for face in self.patches.values():
            face.post_process_face()

        self.reduction_rate = self.get_reduction_rate()
        print ("self.reduction_rate")
        print (self.reduction_rate)

    def get_patches(self):
        """
        this function creates a dictionary with the patches objects
        """

        patches = {}

        patches_list = self.parse_boundary_phi_files()

        for patch in patches_list:
            patches [patch['face_id']] = PatchClass(patch,self.path)


        return patches

    def get_reduction_rate(self):
        """
        this function calculates the reduction rate of the isotope in the cfd
        it sums all the outlet Td sum values - therefore multiple inlets or
        outlets are not supported yet
        """
        atom_in = 0
        atom_out = 0
        red_rate = 0

        for patch in self.patches.values():
            if patch.face_type == 'inlet':
                atom_in += abs(patch.post_process_td_flow[-1])
            elif patch.face_type == 'outlet':
                atom_out += abs(patch.post_process_td_flow[-1])


        if atom_in != 0:
            red_rate = atom_out/atom_in

        return red_rate




    def parse_boundary_phi_files(self):
        """
        this function examines phi files and returns a list of dictionary with
        the face_name and type
        at the moment it is not optimized for speed
        """

        #print ("parsing faces...")

        faces = []

        poly_mesh_folder = os.path.join(self.path,'constant','polyMesh')
        boundary_file = os.path.join(poly_mesh_folder,'boundary')

        try:
            inp_file = open(boundary_file,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open boundary file")
        with inp_file:

            text = inp_file.read()

            face_def_pat = re.compile(r"\d+[\n\r\s]+?\(.*?[\n\r\s]+?\)",
                                    re.MULTILINE | re.DOTALL )
            face_number_pat = re.compile(r"(\d+)[\n\r\s]+?\(.*?\)",
                                       re.MULTILINE | re.DOTALL )
            boundary_pat = re.compile(r"[^\s]+[\n\r\s]+?\{.*?\}",
                                     re.MULTILINE | re.DOTALL )
            boundary_name_pat = re.compile(r"([^\s]+)[\n\r\s]+?\{.*?\}",
                                         re.MULTILINE | re.DOTALL )
            face_n_pat = re.compile(r"nFaces.*?(\d+)")
            first_face_pat = re.compile(r"startFace.*?(\d+)")
            def_block = face_def_pat.findall(text)[0]
            face_number = int(face_number_pat.findall(text)[0])
            boundary_blocks = boundary_pat.findall(def_block)

            for block in boundary_blocks:
                b_dict = {}
                b_dict['face_id'] = boundary_name_pat.findall(block)[0]
                b_dict['face_elements_n'] = int(face_n_pat.findall(block)[0])
                b_dict['face_first_element'] = int(first_face_pat.findall(block)[0])
                faces.append(b_dict)

            if len(faces) != face_number:
                raise ValueError('Error with the number of faces')



        phi_file_path = os.path.join(self.path,str(self.last_time),'phi')
        try:
            inp_file = open(phi_file_path,'r',encoding="utf8", errors='ignore')
        except IOError:
            print("couldn't open phi file")
        with inp_file:

            face_phi_pat = re.compile(r"\((.{1,}?)\)",
                    re.MULTILINE | re.DOTALL )
            text = inp_file.read()
            wall_face_pat = re.compile(r"value\s+uniform\s+0")

            #print ("face phi")
            for face in faces:
                face_block_pat=re.compile(face['face_id'] + r"[\n\r\s]+?\{.*?\}",
                                         re.MULTILINE | re.DOTALL )
                face_block = face_block_pat.findall(text)[0]
                #print (faceBloc)
                face_phis = face_phi_pat.findall(face_block)
                #print (facePhis)
                wall_confirm = wall_face_pat.findall(face_block)
                #print (wallConfirm)
                if len(face_phis) == 0 and len(wall_confirm)!=0 :
                    face['type'] = "wall"
                    face['sum_phis'] = 0
                else:
                    phi_list=face_phi_pat.findall(face_block)[0].strip().split('\n')
                    phi_list = np.array([float(val) for val in phi_list])


                    if all(i >= 0 for i in phi_list):
                        face['type'] = 'outlet'
                        face['sum_phis'] = sum(phi_list)
                    elif all(i <= 0 for i in phi_list):
                        face['type'] = 'inlet'
                        face['sum_phis'] = sum(phi_list)
                    else:
                        raise ValueError('Error, phis with mixed sign in boundary')


        #for ls in faces:
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


        if dimension_vec == [ 1, 0, -1, 0, 0, 0, 0]:
            volumetric_flag = False
        elif dimension_vec == [ 0, 3, -1, 0, 0, 0, 0]:
            volumetric_flag = True
        else:
            raise ValueError('The dimensions of phi are not recognized')

        return volumetric_flag

    def get_density(self):
        """
        this function gets the density value from the
        constant/transportProperties file
        """

        file_path = self.path + '/constant/transportProperties'
        dict_path = 'rhoRef'

        value = self.query_of_single_value(file_path, dict_path)


        return float(value)

    def get_phi_dimensions(self):
        """
        this function gets the dimensions of the phi values from the last
        time folder
        """

        file_path = self.path + '/' + str(self.last_time) + '/phi'

        dimension_vec = self.query_dimensions(file_path)

        return dimension_vec

    def query_dimensions(self, file_path):
        """
        this function queries the phi file to get the dimensions of the phi
        values
        """

        dict_path = 'dimensions'

        dim_vector = []

        for tok_dict in self.token_classifier(file_path):
            if tok_dict['path'] == dict_path and tok_dict['type'] == 'VECTOR_DATA':
                dim_vector.append(int(tok_dict['token']))
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
            if tok_dict['path'] == dict_path:
                field_vector.append(tok_dict['token'])
                found = True

            if found and tok_dict['path'] != dict_path:
                break

        #print (field_vector)

        return field_vector[position]

    def tokenizer(self,file_path):
        """
        generator of the OF file tokens
        """
        i = 1

        with open(file_path, "r",encoding='utf-8') as file:
            for line in file:
                line = line.strip("\r\n")
                line_tokens = re.finditer(
                    r"[^\s{};()/*<>\[\]]+|{|}|;|\(|\)|//|/\*|\*/|\*|<|>|\[|\]|/|\s+", line
                )
                # print (line_tokens)

                #if i > 40:
                #    break
                for token in line_tokens:
                    yield (i, token.group())

                i += 1

    def token_classifier(self,file_path):
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

        for line, token in self.tokenizer(file_path):
            #print("line: ", line, "current_line: ", current_line, "token: ", token)
            j += 1
            token_dict = {}
            # new line
            if line != current_line:
                #output_file.write(line_string + "\n")
                #line_string = ""
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
                #print (token_dict)
                #line_string += token
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

        with open(out_path, "w", encoding='utf-8') as output_file:
            line_string = ""
            current_line = 1
            for token_dict in token_generator:
                # new line
                if token_dict['line'] != current_line:
                    output_file.write(line_string + "\n")
                    line_string = ""
                    current_line = token_dict['line']
                line_string += token_dict['token']

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
        i_density = 1/self.density_kg_m3

        file_path = self.path + '/' + str(self.last_time) + '/phi'

        for t_dict in self.token_classifier(file_path):

            # convert dimensions vector
            if t_dict['path'] == 'dimensions' and t_dict['type'] == 'VECTOR_DATA' :
                count_vec += 1
                if count_vec == 1:
                    t_dict['token'] = str(0)
                if count_vec == 2:
                    t_dict['token'] = str(3)

            if t_dict['path'] == 'internalField' and t_dict['type'] == 'NUMEROSITY':
                internal_check = int(t_dict['token'])
                found_internal = True

            if t_dict['path'] == 'internalField' and t_dict['type'] == 'LIST_DATA':
                count_internal += 1
                new_phi = float(t_dict['token'])*i_density
                t_dict['token'] = f'{new_phi:.11e}'

            if found_internal and t_dict['type'] == 'ROUND_CLOSE':
                if internal_check != count_internal:
                    raise ValueError('ERROR with the conversion of the \
                                     internal phi values')
                found_internal = False

            if ('boundaryField' in t_dict['path'] and
                'value' in t_dict['path'] and
                 t_dict['path'] != name_boundary) :
                found_boundary = True
                name_boundary = t_dict['path']

            if found_boundary and t_dict['token'] == 'uniform':
                single_val = True


            if found_boundary and t_dict['type'] == 'NUMEROSITY':
                boundary_check = int(t_dict['token'])

            if found_boundary and t_dict['type'] == 'LIST_DATA':
                count_boundary += 1
                new_phi = float(t_dict['token'])*i_density
                t_dict['token'] = f'{new_phi:.11e}'

            if found_boundary and t_dict['type'] == 'ROUND_CLOSE':
                if boundary_check != count_boundary:
                    raise ValueError('ERROR with the conversion of the \
                                      boundary phi values')
                count_boundary = 0

            if found_boundary and single_val and is_float(t_dict['token']):
                if t_dict['token'] != '0':
                    new_phi = float(t_dict['token'])*i_density
                    t_dict['token'] = f'{new_phi:.11e}'
                single_val = False

            if found_boundary and t_dict['type'] == 'CURLY_CLOSE':
                found_boundary = False

            yield t_dict