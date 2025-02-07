"""
this file containes the class of the circuit object
"""
#import math
import os
import sys
import numpy as np
from .circuit_node import CircuitNode
from .utils import get_circuit_files, split_text_file_lines,is_float


class CircuitObject:
    """
    this class contains the circuit object
    """
    def __init__(self, item, parameters):
        """
        initialize circuit object
        """

        self.name = item[0]
        self.full_path = os.path.abspath(item[1])
        self.circuit_files = get_circuit_files(self.full_path)

        # create a dictionary that contains the node objects
        # by parsing the nodes file
        self.circuit_dictionary = self.parse_nodes_file()



        # update the circuit dictionary with the info in links.dat
        self.parse_links_file()


        # update the circuit dictionary with the info in the reaction rate mesh (if present)
        self.probe_rrmesh(parameters)

        # create a dictionary that contains the external sources
        # by parsing the the inlets file. the external sources are defined
        # with a simple dictionary
        self.external_nodes = self.define_external_sources()




        # for key,value in self.circuitDictionary.items():
        #    value.printNode()




    def parse_nodes_file(self):
        """
        this function update the circuit dictionary of the circuitObject
        class with the info in the nodes.dat file
        """

        node_dictionary = {}

        nodes_lsplit = split_text_file_lines(self.circuit_files["nodesPath"])

        for line in nodes_lsplit:
            if not line[0].isdigit():
                continue

            node_parameters = {}

            node_parameters['axis_x']      = 0
            node_parameters['axis_y']      = 0
            node_parameters['axis_z']      = 0
            node_parameters['origin_x_cm'] = 0
            node_parameters['origin_y_cm'] = 0
            node_parameters['origin_z_cm'] = 0
            node_parameters['length_cm']   = 0
            node_parameters['radius_cm']   = 0
            node_parameters['axis_x_1'] = 0
            node_parameters['axis_y_1'] = 0
            node_parameters['axis_z_1'] = 0

            node_parameters['node_id'] = int(line[0])
            node_parameters['group_name'] = line[1]
            node_parameters['node_circuit'] = self.name
            node_parameters['node_type'] = line[2].lower()

            node_parameters['node_rtd_time'] = []
            node_parameters['node_rtd_data'] = []
            node_parameters['node_cfd_path'] = ''

            node_parameters['node_stl_path'] = ''

            # determination of the node type value
            # in general it can be a number,
            # a path to a rtd distribution
            # or a path to a fluned simulation
            if is_float(line[3]):
                node_parameters['node_type_value'] = float(line[3])

            elif node_parameters['node_type'] not in ['cfd']:
                # node type value parameter (time treatme ) set to -2 when a
                # residence time distribution is provided
                node_parameters['node_type_value'] = -2
                data_path=os.path.join(self.full_path,line[3].strip('\"\''))
                rtd_node_split = split_text_file_lines(data_path)
                time_vec = []
                data_vec = []
                for  rtd_line in rtd_node_split:
                    time_vec.append(float(rtd_line[0]))
                    data_vec.append(float(rtd_line[1]))
                node_parameters['node_rtd_time'] = np.array(time_vec)
                node_parameters['node_rtd_data'] = np.array(data_vec)
            elif node_parameters['node_type'] in ['cfd']:
                node_parameters['node_type_value'] = -2
                data_path=os.path.join(self.full_path,line[3].strip('\"\''))
                node_parameters['node_cfd_path'] = data_path

            # if the node is stl type look for the stl file path
            if node_parameters['node_type'].lower() == 'stl':
                # find files in the folders that ends with the node.stl extension
                stl_files = []
                for root, _, files in os.walk(self.full_path):
                    for file in files:
                        if file.lower().endswith(f"{node_parameters['node_id']}.stl"):
                            stl_files.append(os.path.join(root, file))
                if len(stl_files) != 1:
                    print(f"ERROR: found {len(stl_files)} stl files for node {node_parameters['node_id']}")
                else:
                    node_parameters['node_stl_path'] = os.path.join(root, stl_files[0])




            # assign physical parameters
            node_parameters['activity_scaling'] = float(line[4])
            node_parameters['temperature_k'] = float(line[5])
            node_parameters['pressure_mpa'] = float(line[6])
            node_parameters['mcnp_cell'] = int(line[7])
            node_parameters['volume_cm3'] = float(line[8])
            node_parameters['reaction_rate_cm3'] = float(line[9])



            # assign cylinder geometry parameters
            if node_parameters['node_type'] in ['pipe', 'tank-cyl']:

                node_parameters['axis_x'] = float(line[10])
                node_parameters['axis_y'] = float(line[11])
                node_parameters['axis_z'] = float(line[12])
                node_parameters['origin_x_cm'] = float(line[13])
                node_parameters['origin_y_cm'] = float(line[14])
                node_parameters['origin_z_cm'] = float(line[15])
                node_parameters['length_cm'] = float(line[16])
                node_parameters['radius_cm']= float(line[17])

                # define axis for theta = 0
                vec_1 = np.array([
                                  float(line[10]),
                                  float(line[11]),
                                  float(line[12]),
                                 ])
                vec_2 = np.random.rand(3)
                vec_2 = vec_2 - np.dot(vec_2, vec_1)/np.dot(vec_1, vec_1)*vec_1
                vec_2 = vec_2/np.linalg.norm(vec_2)
                node_parameters['axis_x_1'] = vec_2[0]
                node_parameters['axis_y_1'] = vec_2[1]
                node_parameters['axis_z_1'] = vec_2[2]


            new_node = CircuitNode(**node_parameters)

            node_dictionary[node_parameters['node_id']] = new_node

        return node_dictionary

    def probe_rrmesh(self, parameters):
        """
        this function updates the circuit dictionary with the info in the
        obtained by probing the rr mesh file
        """

        print ("debug: circuit files", self.circuit_files)

        if "rrmeshPath" not in self.circuit_files:
            return

        for node in self.circuit_dictionary.values():
            node.probe_rrmesh(self.circuit_files["rrmeshPath"], parameters)

        return 0


    def parse_links_file(self):
        """
        this function parses the links file and add the parents to the node
        objects
        """

        link_lsplit = split_text_file_lines(self.circuit_files["linksPath"])

        for i, line in enumerate(link_lsplit):
            # print (line)
            if not line[0].isdigit():
                continue

            node_id = int(line[0])
            #n_parents = int(line[1])
            circuit_name = self.name

            parents_vec = []

            if len(line) == 2:
                # parent definition is not present: the parent is the cell
                # in the line above
                parent_id = 1
                parent_vec = []
                parent_vec.append(parent_id)
                parent_vec.append("int")
                default_parent_id = link_lsplit[i - 1][0]
                parent_vec.append(default_parent_id)
                parent_vec.append(1.00)
                #print (parent_vec)

                parents_vec.append(parent_vec)

            if len(line) > 2:
                # parent definition is in the line
                parent_id = 1
                parent_vec = []
                parent_vec.append(parent_id)
                for item in line[2:]:
                    parent_vec.append(item)
                    if len(parent_vec) == 4:
                        #print (parent_vec)
                        parents_vec.append(parent_vec)
                        parent_id += 1
                        parent_vec = []
                        parent_vec.append(parent_id)

            for parent in parents_vec:
                parent_node_id = int(parent[2])
                parent_type = parent[1]
                child_fraction = float(parent[3])
                # if the parent is an internal node we can already add the
                # parent and child information to the corresponding nodes
                if parent_type == "int":
                    self.circuit_dictionary[node_id].add_parent(parent,circuit_name)
                    self.circuit_dictionary[parent_node_id].add_child(
                                            node_id,
                                            parent_type,
                                            child_fraction,
                                            circuit_name)
                # if the parent is external we can set only the parent information
                # while the child information will be set when all the circuits
                # have been initialized
                elif parent_type == "ext":
                    self.circuit_dictionary[node_id].add_parent(parent)


        return 0

    def define_external_sources(self):
        """
        this function creates an external sources dictionary by parsing
        the inlets.dat file. It also checks the correspondance with the
        circuit dictionary
        """

        inlet_lsplit = split_text_file_lines(self.circuit_files["inletsPath"])

        external_dictionary = {}
        check_inlets = []

        for  line in inlet_lsplit:
            if not line[0].isdigit():
                continue

            inlet_dictionary = {}
            in_id = int(line[0])
            inlet_dictionary["ID"] = in_id
            inlet_dictionary["type"] = line[1]
            inlet_dictionary["value"] = line[2]

            if inlet_dictionary["type"] == "circuit-node":
                connection_def = inlet_dictionary["value"].split("/")
                inlet_dictionary["circuit_source"] = connection_def[0]
                inlet_dictionary["circuit_source_node"] = int(connection_def[1])

                # in this case the inlet closes a loop - the last value after
                # the symbol / is the mass inflow
                if len(connection_def) == 3:
                    inlet_dictionary["loop_mass"] = float(connection_def[2])

            if inlet_dictionary["type"] == "source-node":
                inlet_dictionary["mass_inflow"] = float(inlet_dictionary["value"])
                inlet_dictionary["isotope"] = line[3]
                if is_float(line[4]):
                    inlet_dictionary["isotope_value"] = float(line[4])
                else:
                    data_path=os.path.join(self.full_path,line[4].strip('\"\''))
                    inlet_act_split = split_text_file_lines(data_path)
                    time_vec = []
                    data_vec = []
                    for  line in inlet_act_split:
                        time_vec.append(float(line[0]))
                        data_vec.append(float(line[1]))

                    inlet_dictionary["isotope_value"] = -2
                    inlet_dictionary["isotope_value_time"] = time_vec
                    inlet_dictionary["isotope_value_data"] = data_vec


                if inlet_dictionary["isotope"] not in ["N16", "N17", "O19"]:
                    raise ValueError("ERROR: isotope not implemented - EXITING")

            check_inlets.append(in_id)
            external_dictionary[in_id] = inlet_dictionary


        # check the correspondance between the two input files
        # get the external IDs from the circuit Dictionary
        check_dict = []
        for item in self.circuit_dictionary.values():
            for parent in item.parents.values():
                if parent["parent_type"] == "ext":
                    check_dict.append(parent["parent_node_id"])

        check_dict = sorted(list(set(check_dict)))
        check_inlets = sorted(list(set(check_inlets)))


        if check_dict != check_inlets:
            print("ERROR the external sources defined in:")
            print("inlets.dat and links.dat do not match! EXITING")
            raise ValueError()

        return external_dictionary
