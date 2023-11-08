"""
class for the OF simulations, this can be used to parse and generate files
"""
import os
import re
import sys

def is_float(s):
    """
    check if a string can be converted to float
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

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