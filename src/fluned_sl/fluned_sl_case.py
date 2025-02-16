"""
this file containes the class of the simulation object
"""
import os
import sys
import math

import numpy as np
from scipy.sparse import csc_matrix

from scipy.sparse.linalg import spsolve

from .circuit_object import CircuitObject
from .utils import unpack_mcnp, search_circuits, interp_series

from .pipes_cdgs import make_pipes_cdgs

from .pipes_vtk import write_vtk_files

from water_isotopes.water_isotopes import get_isotope_data

from .tdd_functions import (build_augmented_matrix,
                            build_forcing_vector,
                            build_rr_augmented_matrix,
                            build_decay_average_augmented_matrix,
                            build_rr_average_augmented_matrix)

class FlunedSlCase:
    """
    this class contains the main case object
    """

    def __init__(self, user_parameters, case_path):
        """
        initialize case, look for circuits input and
        initialize these
        """

        default_parameters = {
            "steady_state": True,
            "isotope" : "n16",
            "fluid_carrier": "water",
            "n_sampling_min": 1e05,
            "n_sampling_max": 1e08,
            "rrmesh_sampling_error_max":0.20,
            "rrmesh_sampling_field_name":"Value - Total",
            "rrmesh_sampling_errfield_name" : "Error - Total",
            "rrmesh_sampling_cm":2.0,
            "rrmesh_scaling_factor":1.0,
            "rrmesh_time_scaling_series": '',
            "pipe_time_default": "reynolds",
            "tank_time_default": "uniform",
            "tank_cyl_time_default": "uniform",
            "t_delta_sec": 1.0,
            "t_delta_sec_plot": 1.0,
            "t_begin_sec": 0.0,
            "t_end_sec":  10.0,
            "mc_error_max": 1e-2,
            "t_bins": "0 1I 2",
            "numerical_method": "monte_carlo",
        }

        for key, value in default_parameters.items():
            if key not in user_parameters:
                user_parameters[key] = value


        self.parameters = user_parameters

        self.case_directory = case_path

        for key, value in default_parameters.items():
            if key not in user_parameters:
                setattr(self, key,value)
            else:
                setattr(self, key,user_parameters[key])



        self.single_isotope = True

        # this variable contains a list of tuples where the first item is the
        # the folder name (which is used as the circuit name) and the second
        # item is the full path to the folder
        self.circuit_directories = search_circuits(self.case_directory)

        # initialize sub-circuits
        self.circuits = {}
        for item in self.circuit_directories:
            self.circuits[item[0]] = CircuitObject(item, self.parameters)


        # update the children and parent dictionaries with the info in inlet.dat
        self.define_circuit_links()

        self.all_nodes = self.get_all_nodes()
        self.tot_n_nodes = len(self.all_nodes)



        if self.single_isotope:

            self.isotope_decay_constant = 0.0
            self.isotope_branching_ratio = 0.0
            self.isotope_particle = ''
            self.isotope_energy_bins = []
            self.isotope_energy_probs  = []

            self.define_single_isotope_properties()


        # calculate mass flows
        self.calculate_mass_flows()

        self.define_parents_weights()

        # now that the child dictionaries are defined for each node - calculate
        # the probability vector for each node
        self.define_children_weights()

        # assign density
        if self.fluid_carrier == "water":
            self.assign_density_water()
        else:
            print ("fluid carrier not recognized")
            sys.exit()

        # calculate bulk velocity and bulk residence time for each node
        self.calculate_bulk_properties()

        # initialize the array of independent counters
        self.ext_source_list,self.ext_source_prob, self.mc_list, self.source_idx = self.define_start_nodes()
        self.n_sources = len(self.mc_list)
        self.reset_node_counters()

        # load the reaction rate scaling series
        self.rr_scaling_series = self.load_rr_scaling_series()

    def load_rr_scaling_series(self):
        """
        this function loads the reaction rate scaling series
        """

        if self.rrmesh_time_scaling_series == '':
            return []

        rr_scaling_series = {}
        time = []
        values = []


        # try to open the file
        try:
            with open(self.rrmesh_time_scaling_series, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line == '' or not line[0].isdigit():
                        continue
                    line = line.split(',')
                    time.append(float(line[0]))
                    values.append(float(line[1]))
        except:
            print ("ERROR: could not open the rr scaling series file")
            sys.exit()

        rr_scaling_series['time'] = time
        rr_scaling_series['values'] = values




        return rr_scaling_series


    def define_single_isotope_properties(self):
        """
        In case of a single isotope run this function get the isotope properties
        """


        if self.isotope == "N17":
            self.isotope_particle = "neutron"
        else:
            self.isotope_particle = "photon"

        isotope_database = get_isotope_data()

        isotope_data = isotope_database[self.isotope]


        self.isotope_decay_constant    = isotope_data['decay_constant']
        self.isotope_branching_ratio   = isotope_data['branching_ratio']
        self.isotope_energy_bins     = isotope_data['e_bins']
        self.isotope_energy_probs    = isotope_data['p_bins']

        return 0


    def assign_density_water(self):
        """
        this function assigns the density to the nodes
        assuming it is water. It uses the iapws model based on the
        IAPWS-IF97 formulation
        """

        for node_keys in self.all_nodes:
            node = self.get_node(node_keys)
            node.assign_node_water_properties()

        return 0

    def calculate_bulk_properties(self):
        """
        once a density and a mass flow is assigned to each node we can
        calculate its bulk properties:

        reynolds
        bulk velocity
        uniform velocity - residence time
        uniform velocity - decay rate
        uniform velocity - normalized average activity

        """

        for node_keys in self.all_nodes:
            node = self.get_node(node_keys)
            node.calculate_bulk_properties()
            node.calculate_bulk_activation(self.isotope_decay_constant)
            node.preprocess_rtd_data(self.isotope_decay_constant)

        return 0

    def define_circuit_links(self):
        """
        this function takes the definition of the source nodes and add the
        corresponding information to the node object dictionary for the link
        nodes
        """

        # check links between circuits consistency
        self.check_links()

        for circuit in self.circuits.values():
            for ext_source in circuit.external_nodes.values():


                if ext_source['type'] == 'source-node':

                    # update the parent dictionary of the child node
                    for node in circuit.circuit_dictionary.values():
                        for parent in node.parents.values():

                            if (parent['parent_type'] == 'ext' and
                                parent['parent_node_id'] == ext_source['ID']):

                                parent['parent_ext_type'] = 'source'
                                parent['parent_ext_mass_inflow'] = (
                                            ext_source['mass_inflow'])
                                parent['parent_ext_activity_inflow'] = (
                                            ext_source['isotope_value'])
                                if parent['parent_ext_activity_inflow'] == -2 :
                                    parent['parent_ext_activity_data'] = (
                                            ext_source['isotope_value_data'])
                                    parent['parent_ext_activity_time'] = (
                                            ext_source['isotope_value_time'])
                                ext_source['child_keys'] = ((circuit.name,
                                                            node.node_id))
                                break

                if ext_source['type'] == 'circuit-node':

                    # update the parent dictionary of the child node
                    for node in circuit.circuit_dictionary.values():
                        for parent in node.parents.values():
                            if (parent['parent_type'] == 'ext' and
                                parent['parent_node_id'] == ext_source['ID']):

                                child_fraction = parent['fraction']
                                parent_circuit = ext_source['circuit_source']
                                parent_node_id = ext_source['circuit_source_node']
                                ext_type = 'node'

                                if "loop_mass" in ext_source:
                                    parent['parent_ext_mass_inflow'] = (
                                            ext_source['loop_mass'])
                                    ext_type = 'loop'


                                parent['parent_ext_type'] = 'node'
                                parent['parent_circuit'] = parent_circuit
                                parent['parent_keys'] = (parent_circuit,
                                                         parent_node_id)

                                child_node_id = node.node_id
                                child_circuit = circuit.name

                                parent_node = self.get_node([parent_circuit,
                                                            parent_node_id])

                                parent_node.add_child(child_node_id,
                                                      'ext',
                                                       child_fraction,
                                                       child_circuit,
                                                       ext_type = ext_type)



        return 0

    def check_links(self):
        """
        this function run a check that the links between circuits
        are valid
        """

        for circuit_value in self.circuits.values():
            for ext_value in circuit_value.external_nodes.values():
                if ext_value["type"] == "circuit-node":
                    test_key_1 = ext_value["circuit_source"]
                    test_key_2 = ext_value["circuit_source_node"]
                    test_link = self.get_node((test_key_1, test_key_2))
                    if not test_link:
                        print("ERROR link not found:")
                        print("circuit: ", test_key_1)
                        print("node: ", test_key_2)

        return 0

    def define_children_weights(self):
        """
        Once all the node links are defined this function calculates the
        probability vector for each branch
        """

        for circuit in self.circuits.values():

            for node in circuit.circuit_dictionary.values():

                node.define_children_weights()

        return 0

    def define_parents_weights(self):
        """
        Once all the node links are defined this function calculates the actual
        parent weights for each node
        """
        for circuit in self.circuits.values():
            for node in circuit.circuit_dictionary.values():
                frac_sum = 0
                node_keys = node.node_keys
                for parent in node.parents.values():

                    # connection to internal line
                    if parent['parent_type'] == 'int':
                        parent_keys = parent['parent_keys']
                        parent_node = self.get_node(parent_keys)
                        parent_inflow = parent_node.mass_flow*parent['fraction']
                        parent['parent_int_mass_inflow'] = parent_inflow
                        parent['fraction'] = parent_inflow/node.mass_flow

                        for child in parent_node.children.values():
                            if child['child_keys'] == node_keys:
                                child['fraction_dstream'] = parent['fraction']
                                break


                    elif parent['parent_type'] == 'ext' and parent['parent_ext_type'] == 'node':
                        parent_keys = parent['parent_keys']
                        parent_node = self.get_node(parent_keys)
                        parent_inflow = parent_node.mass_flow*parent['fraction']
                        parent['parent_ext_node_mass_inflow'] = parent_inflow
                        parent['fraction'] = parent_inflow/node.mass_flow
                        for child in parent_node.children.values():
                            if child['child_keys'] == node_keys:
                                child['fraction_dstream'] = parent['fraction']
                                break
                    elif parent['parent_type'] == 'ext' and parent['parent_ext_type'] == 'source':
                        parent['fraction'] = parent['parent_ext_mass_inflow']/node.mass_flow

                    frac_sum += parent['fraction']

                if not math.isclose(frac_sum, 1, rel_tol=1e-5):
                    print (f"ERROR: parent fractions of node {node.node_id} in group {node.group_name} do not sum to 1")
                    print (f"frac sum value: ", frac_sum)
                    sys.exit()

        return 0


    def loop_check(self):
        """
        this function checks if there are loops in the circuit. Meaning that
        at some point it goes back to the same node. At the moment this is not
        implemented for steady state simulations.
        """
        print ("circuit sampling for loops")

        i = 0

        step = pow(10,int(math.log10(self.n_sampling_min)-1))


        check = False
        loop_check = False

        while True:


            # get a random start node
            start_node, _ ,_= self.pick_start_node()

            # traverse the whole circuit until the end node
            traverse_nodes, _ = self.traverse_nodes_forward(start_node)

            if len(traverse_nodes) < len(set(traverse_nodes)):
                print("loop detected")
                print(traverse_nodes)
                loop_check = True
                break

            for node_keys in traverse_nodes:
                node = self.get_node(node_keys)
                node.increase_mc_counter()

            i += 1

            if (i) % step == 0:
                check = self.eval_mc_checkpoint(self.mc_error_max,
                                                self.n_sampling_min,i)
                if check:
                    break

            if i >= self.n_sampling_max:
                break


        self.reset_node_counters()

        if loop_check:
            print ("loop detected in the circuit")
        else:
            print ("no loop detected in the circuit")


        return loop_check

    def calculate_atom_flows(self):
        """
        define and solve a balance matrix and later calculates the average
        and outlet activity
        """

        matrix_list = self.all_nodes

        m_size = len(matrix_list)
        b_vector = np.zeros(m_size)

        data = []
        row_ind = []
        col_ind = []

        for i, elem in enumerate(matrix_list):

            row_node = self.get_node(elem)
            #density_kg_m3 = row_node.density_g_cm3*1000
            #b_mass = row_node.get_external_inflow()
            b_act = row_node.get_external_inflow_activity()
            b_vector[i] += b_act

            # append diagonal elements
            data.append(1)
            row_ind.append(i)
            col_ind.append(i)


            for child in self.get_node(elem).children.values():
                child_keys = child["child_keys"]
                child_mflow = self.get_node(child_keys).mass_flow
                j = matrix_list.index(child_keys)
                parent_mflow = row_node.mass_flow
                parent_mflow_fraction = child["fraction"]

                out_activity_factor, _ = row_node.get_outlet_average_activity(
                                            1,
                                            self.tank_time_default,
                                            self.tank_cyl_time_default,
                                            self.pipe_time_default,
                                            self.isotope_decay_constant,
                                            steady_state = True,
                                            activation_flag = False,
                                            numerical_method = "deterministic",
                )

                act_flow_dstream = -(parent_mflow*
                                    parent_mflow_fraction*
                                    out_activity_factor/
                                    child_mflow)

                #if i  in [0,1]:
                #    print ("parent: ", elem)
                #    print("child: ", child_keys)
                #    print ("mass flow: ", parent_mflow)
                #    print ("mass flow fraction: ", parent_mflow_fraction)
                data.append(act_flow_dstream)
                row_ind.append(j)
                col_ind.append(i)

                out_activity_rr, _ = row_node.get_outlet_average_activity(
                                            0,
                                            self.tank_time_default,
                                            self.tank_cyl_time_default,
                                            self.pipe_time_default,
                                            self.isotope_decay_constant,
                                            steady_state = True,
                                            activation_flag = True,
                                            numerical_method = "deterministic",
                )
                rr_flow = (parent_mflow*
                          parent_mflow_fraction*
                          out_activity_rr/
                          child_mflow)

                #if i  in [0,1]:
                #    print ("rr flow: ", rr_flow)
                #    print("child rr: ", j)
                b_vector[j] += rr_flow




        a_flow_matrix = csc_matrix((data,
                                   (row_ind, col_ind)),
                                   shape=(m_size, m_size))

        #print (a_flow_matrix.toarray())
        #print ("b vector: ", b_vector)

        aflow = spsolve(a_flow_matrix, b_vector)

        #print ("a flow: ", aflow)

        rtol = 1e-5
        rel_diff = ([abs(a-b)/b if b != 0 else 0 for a, b in
                     zip(a_flow_matrix.dot(aflow), b_vector)])




        if not all(x < rtol for x in rel_diff ):
            print ("a flow matrix: ", a_flow_matrix.toarray())
            print ("b vector: ", b_vector)
            print("ERROR! atom flow balance solution not solved")
            print ("relative difference: ", rel_diff)
            print ("max relative difference: ", max(rel_diff))
            sys.exit()


        # assign activity flow to nodes
        for i, keys in enumerate(matrix_list):
            node = self.get_node(keys)
            #density_kg_m3 = node.density_g_cm3*1000
            act_in = aflow[i]#*density_kg_m3/node.mass_flow
            #print ("node: ", node.node_keys)
            #print ("act in solution: ", aflow[i])
            #print ("node density: ", density_kg_m3)
            #print ("node mass flow: ", node.mass_flow)
            #print ("act in: ", act_in)
            #print ("mass flow: ", node.mass_flow)
            #print ("residence time: ", node.sample_res_time)
            #print ("source term: ", b_vector[i])
            act_out, act_average = node.get_outlet_average_activity(
                                            act_in,
                                            self.tank_time_default,
                                            self.tank_cyl_time_default,
                                            self.pipe_time_default,
                                            self.isotope_decay_constant,
                                            steady_state = True,
                                            activation_flag = True,
                                            numerical_method = "deterministic",
                                        )
            #print ("act out: ", act_out)
            #print()
            node.inlet_activity_bq_m3 = act_in
            node.out_activity_bq_m3 = act_out
            node.average_activity_bq_m3 = act_average
            vol_m3 = node.volume_cm3 * 1e-6
            node.tot_activity_bq = act_average*vol_m3

        return 0


    def calculate_mass_flows(self):
        """
        define and solve a balance mass matrix and later assign all the
        mass flows to the nodes in kg/s
        """

        matrix_list = self.all_nodes

        m_size = len(matrix_list)
        b_vector = np.zeros(m_size)

        data = []
        row_ind = []
        col_ind = []

        for i, elem in enumerate(matrix_list):

            row_node = self.get_node(elem)
            sum_b = row_node.get_external_inflow()
            b_vector[i] = -sum_b

            # append diagonal elements
            data.append(-1)
            row_ind.append(i)
            col_ind.append(i)


            for child in self.get_node(elem).children.values():
                child_keys = child["child_keys"]
                j = matrix_list.index(child_keys)
                if child["child_ext_type"] == "loop":
                    continue
                data.append(child["fraction"])
                row_ind.append(j)
                col_ind.append(i)

        m_flow_matrix = csc_matrix((data,
                                   (row_ind, col_ind)),
                                   shape=(m_size, m_size))

        #print (b_vector)
        #print (m_flow_matrix.toarray())
        mflow = spsolve(m_flow_matrix, b_vector)

        # print (exitCode)
        if not np.allclose(m_flow_matrix.dot(mflow), b_vector):
            print ("m flow matrix: ", m_flow_matrix.toarray())
            print ("b vector: ", b_vector)
            print("ERROR! mass flow balance solution not solved")
            sys.exit()

        # assign mass flow to nodes
        for i, keys in enumerate(matrix_list):
            node = self.get_node(keys)
            node.mass_flow = mflow[i]

        return 0

    def solve(self):
        """
        this function calculates the fluid activation
        """

        if self.single_isotope:

            if self.steady_state:

                if self.numerical_method == "monte_carlo":
                    self.mc_solve_single_isotope_steady_state()
                elif self.numerical_method == "deterministic":
                    self.det_solve_single_isotope_steady_state()

            else:

                if self.numerical_method == "monte_carlo":
                    self.mc_solve_single_isotope_transient()
                elif self.numerical_method == "deterministic":
                    self.det_solve_single_isotope_transient()

        else:

            print("multiple isotopes not implemented yet")
            sys.exit()

        return 0


    def mc_solve_single_isotope_transient(self):
        """
        this function solves the monte carlo problem for a single isotope
        in transient mode
        -
        activation inside the circuit is not present at the moment
        """

        print("solving transient single isotope problem")


        mc_step = pow(10,int(math.log10(self.n_sampling_min)-1))
        time_steps = unpack_mcnp(self.parameters["t_bins"])

        for t_step in time_steps:

            if t_step == 0:
                self.write_vtk_single_isotope(0)
                continue

            self.reset_node_counters()

            print ("time step: ", t_step)

            j = 0
            check = False

            while True:


                # get a random start node
                start_node, _, _ = self.pick_start_node()

                #calculate a random radius
                radius_samp = np.sqrt(np.random.uniform(0.0,1.0))

                # traverse the whole circuit until the end node
                traverse_nodes, _ = self.traverse_nodes_forward(start_node)
                mc_id = self.pick_start_mc_counter()


                # calculate a vector with the cumulated residence time along
                # the sampled path. In this way it is possible to cut the
                # calculation path where the wavefront is


                t_res_vec = np.zeros(len(traverse_nodes))
                t_res_complete_vec = np.zeros(len(traverse_nodes))
                t_in_vec = np.zeros(len(traverse_nodes))

                print ("traverse nodes", traverse_nodes)
                print ("mc id", mc_id)

                res_time_in = 0.0
                found_source = False
                for i,node_keys in enumerate(traverse_nodes):
                    node = self.get_node(node_keys)
                    found_source = node.set_residence_time(
                                              radius_samp,
                                              self.isotope_decay_constant,
                                              self.tank_time_default,
                                              self.tank_cyl_time_default,
                                              self.pipe_time_default,
                                              found_source_flag=found_source,
                                              mc_counter = mc_id,
                                              )
                    t_res_vec[i] = node.sample_res_time
                    t_res_complete_vec[i] = node.sample_res_time/node.sample_a_frac
                    t_in_vec[i] = res_time_in
                    res_time_in += node.sample_res_time
                    if res_time_in > t_step:
                        traverse_nodes = traverse_nodes[:i+1]
                        break

                print ("traverse nodes cut", traverse_nodes)
                print ("residence time", t_res_vec)
                print ("t in vec", t_in_vec)
                print ("t res complete vec", t_res_complete_vec)
                print ()

                #print ("residence time: ", t_res_vec)
                #print ("inlet time: ", t_in_vec)



                for k, (t_in, node_keys) in enumerate(zip(t_in_vec,traverse_nodes)):
                    t_source = t_step - t_in
                    inlet_node = self.get_node(traverse_nodes[0])
                    act_in = inlet_node.get_external_inflow_activity(
                                 t_sample = t_source,
                                 mc_counter = mc_id,
                                  )
                    t_inside = t_step - t_in
                    #print ("inlet activity: ", act_in)

                    # an inner loop is needed to calculate the inlet activity
                    # of the node at the fixed time step. This process is
                    # recursive but it is needed if component with a rtd, or
                    # with activation are defined. It is possible to add
                    # optimization



                    for h, node_keys_inner in enumerate(traverse_nodes[:k]):
                        node_inner = self.get_node(node_keys_inner)
                        res_t_tot_fraction = 1

                        if t_in_vec[h] == 0 and node_inner.sample_a_frac != 1:
                            t_inside_h = t_step - t_in_vec[h]
                            if t_inside_h >= t_res_complete_vec[h]:
                                res_t_tot_fraction = 1
                            else:
                                res_t_tot_fraction = t_inside_h/t_res_complete_vec[h]
                        act_in, _ = node_inner.mc_node_calculation(
                                            act_in,
                                            t_in,
                                            self.tank_time_default,
                                            self.tank_cyl_time_default,
                                            self.pipe_time_default,
                                            self.isotope_decay_constant,
                                            res_t_fraction = 1,
                                            res_t_complete_fraction = res_t_tot_fraction,
                                            steady_state = False,
                                            outer_loop = False,
                                            mc_counter = mc_id,
                                            )

                    #print ("end of inner loop")



                    node = self.get_node(node_keys)




                    print ("calculation of res_t_tot_fraction")
                    print ("k", k)
                    print ("t inside", t_inside)
                    print ("t res complete vec", t_res_complete_vec)
                    print ("t_res_complete_vec[k]", t_res_complete_vec[k])

                    res_t_tot_fraction = 1
                    if t_in == 0 and node.sample_a_frac != 1:
                        if t_inside >= t_res_complete_vec[k]:
                            res_t_tot_fraction = 1
                        else:
                            res_t_tot_fraction = t_inside/t_res_complete_vec[k]


                    if k != len(traverse_nodes)-1:
                        res_t_fraction = 1
                    else:
                        #last node
                        if t_inside > t_res_vec[k]:
                            res_t_fraction = 1
                        else:
                            res_t_fraction = t_inside/t_res_vec[k]

                    _,_ = node.mc_node_calculation(
                                        act_in,
                                        t_in,
                                        self.tank_time_default,
                                        self.tank_cyl_time_default,
                                        self.pipe_time_default,
                                        self.isotope_decay_constant,
                                        res_t_fraction = res_t_fraction,
                                        res_t_complete_fraction = res_t_tot_fraction,
                                        steady_state = False,
                                        outer_loop = True,
                                        mc_counter = mc_id,
                                        )


                j += 1
                #print("outer loop: ", j)
                #print()


                if (j) % mc_step == 0:
                    check = self.eval_mc_checkpoint(self.mc_error_max,
                                                    self.n_sampling_min,j)
                    if check:
                        break


                if j >= self.n_sampling_max:
                    break


            self.write_cdgs_single_isotope(t_step)
            self.write_vtk_single_isotope(t_step)
            self.print_solution_single_isotope_transient(t_step)




        return 0

    def det_solve_single_isotope_transient(self):
        """
        this function solves the transport of a single isotope with a
        deterministic  method.
        """

        print("solving transient single isotope with deterministic method")


        # at the moment set the pipe_time_default to 'uniform' to ensure we
        # dismiss the radius sampling in the future we might be able to
        # implement it
        self.pipe_time_default = 'uniform'
        self.tank_cyl_time_default = 'uniform'

        test_list_di = []
        test_list_di_err = []



        # first assign a residence time to each node
        for node_keys in self.all_nodes:
            node = self.get_node(node_keys)
            _ = node.set_residence_time(
                                    0,
                                    self.isotope_decay_constant,
                                    self.tank_time_default,
                                    self.tank_cyl_time_default,
                                    self.pipe_time_default,
                                    numerical_method = "deterministic",
                                    )
            node.mc_res_time = node.sample_res_time
            node.calculate_delay_int(self.t_delta_sec)

            test_list_di.append(node.delay_int)
            test_list_di_err.append(node.delay_int_error)


        print ("delay int: ", test_list_di)
        #print ("delay int err: ", test_list_di_err)
        print ("debub max error: ", max(test_list_di_err))
        print ("debug numbers of int delay, ", len(set(test_list_di)))
        print ("debug of int delay set, ", set(test_list_di))

        t_start = self.t_begin_sec
        t_end = self.t_end_sec
        t_delta = self.t_delta_sec

        time_steps = np.linspace(t_start, t_end, num=int((t_end - t_start) / t_delta) + 1)
        n_time_steps = len(time_steps)



        print ("debug time steps: ", time_steps)

        # gather all nodes
        matrix_list = self.all_nodes
        m_size = len(matrix_list)

        node_weights_matrix = np.zeros((m_size,m_size))
        res_time_matrix = np.zeros((m_size,m_size))
        system_decay_matrix = np.zeros((m_size,m_size))
        delay_matrix = np.zeros((m_size, m_size), dtype=int)
        b_vector = np.zeros(m_size)



        # define the system delay matrix
        for i, elem in enumerate(matrix_list):

            row_node = self.get_node(elem)
            #density_kg_m3 = row_node.density_g_cm3*1000
            #b_mass = row_node.get_external_inflow()
            b_act = row_node.get_external_inflow_activity()
            b_vector[i] += b_act



            for child in self.get_node(elem).children.values():
                child_keys = child["child_keys"]
                child_mflow = self.get_node(child_keys).mass_flow
                j = matrix_list.index(child_keys)
                parent_mflow = row_node.mass_flow
                parent_mflow_fraction = child["fraction"]

                out_activity_factor, _ = row_node.get_outlet_average_activity(
                                            1,
                                            self.tank_time_default,
                                            self.tank_cyl_time_default,
                                            self.pipe_time_default,
                                            self.isotope_decay_constant,
                                            steady_state = False,
                                            activation_flag = False,
                                            numerical_method = "deterministic",
                )

                act_flow_dstream = (parent_mflow*
                                    parent_mflow_fraction/
                                    child_mflow)

                node_weights_matrix[j,i] = act_flow_dstream
                system_decay_matrix[j,i] = out_activity_factor
                delay_matrix[j,i] = row_node.delay_int
                res_time_matrix[j,i] = row_node.sample_res_time

        n_max = np.max(delay_matrix)

        # matrix to store the time-vaying solution
        X = np.zeros((n_time_steps, m_size*n_max))

        # matrix to store the time-vaying node average activity
        Av = np.zeros((n_time_steps, m_size*n_max))

        # matrix to store the time varying reaction rates - Initialization
        S = np.zeros((n_time_steps, m_size*n_max))

        # interpolate the rr time series over the n_time_steps
        rr_rates = np.zeros(m_size)
        for i, elem in enumerate(matrix_list):
            row_node = self.get_node(elem)
            rr_rates[i] = row_node.reaction_rate_m3

        rr_mesh_time_series = np.zeros(n_time_steps)
        for i, time in enumerate(time_steps):
            rr_mesh_time_series[i] = interp_series(self.rr_scaling_series['time'],
                                                   self.rr_scaling_series['values'],
                                                   time)

        for k in range(n_time_steps):
            s_row = np.zeros((0))
            k_list = np.arange(k, k-n_max, -1)
            for k_val in k_list:
                if k_val < 0:
                    s_row = np.concatenate((s_row, np.zeros(m_size)))
                else:
                    s_row = np.concatenate((s_row, rr_rates*rr_mesh_time_series[k_val]))
            S[k,:] = s_row






        print ("debug max delay: ", n_max)

        # augmentation matrix for outlet activity from entering activity
        aug_matrix = build_augmented_matrix(node_weights_matrix,
                                                  delay_matrix,
                                                  res_time_matrix,
                                                  self.isotope_decay_constant)


        # augmentation matrix for outlet activity from rr activity
        rr_aug_matrix = build_rr_augmented_matrix(node_weights_matrix,
                                                  delay_matrix,
                                                  res_time_matrix,
                                                  self.isotope_decay_constant)

        # augmentation matrix for average activity from entering activity
        avg_aug_matrix = build_decay_average_augmented_matrix(
                                                delay_matrix,
                                                res_time_matrix,
                                                self.isotope_decay_constant
        )

        # augmentation matrix for average activity from rr activity
        rr_avg_aug_matrix = build_rr_average_augmented_matrix(
                                                delay_matrix,
                                                res_time_matrix,
                                                self.isotope_decay_constant
        )

        u_vector = build_forcing_vector(b_vector, n_max)

        # define the k=0 solution
        X0 = build_forcing_vector(b_vector, n_max)
        X[0,:] = X0.reshape(-1)



        vtk_id = -1

        for k in range(n_time_steps ):
            print ("time step: ", k, "/", (n_time_steps-1))
            if k!= n_time_steps-1:
                X[k+1, :] = (aug_matrix @ X[k, :].reshape(-1) +
                             rr_aug_matrix @ S[k, :].reshape(-1) +
                            u_vector.reshape(-1))
                Av[k+1, :] = (avg_aug_matrix @ (X[k, :].reshape(-1)) +
                             rr_avg_aug_matrix @ S[k, :].reshape(-1) +
                              u_vector.reshape(-1))

            # assign activity flow to nodes
            for i, keys in enumerate(matrix_list):
                node = self.get_node(keys)

                act_in = X[k,i]
                act_average = Av[k,i]

                # act out is calculated in the old method and will be removed
                act_out, _ = node.get_outlet_average_activity(
                                                act_in,
                                                self.tank_time_default,
                                                self.tank_cyl_time_default,
                                                self.pipe_time_default,
                                                self.isotope_decay_constant,
                                                steady_state = True,
                                                activation_flag = False,
                                                numerical_method = "deterministic",
                                            )
                #print ("act out: ", act_out)
                node.inlet_activity_bq_m3 = act_in
                node.out_activity_bq_m3 = act_out
                node.average_activity_bq_m3 = act_average
                vol_m3 = node.volume_cm3 * 1e-6
                node.tot_activity_bq = act_average*vol_m3

            # output data
            t_step = t_delta*k
            self.print_solution_single_isotope_transient(t_step)

            if t_step // self.t_delta_sec_plot == (vtk_id+1):
                vtk_id += 1
                id_string = f"{vtk_id:d}".zfill(int(math.floor(math.log10(n_time_steps))+1))
                self.write_vtk_single_isotope(id_string)



        return

    def det_solve_single_isotope_steady_state(self):
        """
        this function solves the transport of a single isotope with a
        deterministic  method. i.e. by solving a system of linear equations
        """

        print("solving steady state single isotope with deterministic method")


        # at the moment set the pipe_time_default to 'uniform' to ensure we
        # dismiss the radius sampling in the future we might be able to
        # implement it
        self.pipe_time_default = 'uniform'
        self.tank_cyl_time_default = 'uniform'



        # first assign a residence time to each node
        for node_keys in self.all_nodes:
            node = self.get_node(node_keys)
            _ = node.set_residence_time(
                                    0,
                                    self.isotope_decay_constant,
                                    self.tank_time_default,
                                    self.tank_cyl_time_default,
                                    self.pipe_time_default,
                                    numerical_method = "deterministic",
                                    )
            node.mc_res_time = node.sample_res_time

        # solve the equation system
        self.calculate_atom_flows()





        self.write_vtk_single_isotope()
        self.write_cdgs_single_isotope()
        self.print_solution_single_isotope_steady_state()

        return 0


    def mc_solve_single_isotope_steady_state(self):
        """
        this function solves the monte carlo problem for a single isotope
        in steady state
        """

        print("solving steady state single isotope problem")

        i = 0

        step = pow(10,int(math.log10(self.n_sampling_min)-1))

        check = False

        while True:


            # get a random start node
            start_node, _, _ = self.pick_start_node()

            #sample a random radius
            radius_samp = np.sqrt(np.random.uniform(0.0,1.0))

            # traverse the whole circuit until the end node
            traverse_nodes, _ = self.traverse_nodes_forward(start_node)

            #print ("traverse nodes: ", traverse_nodes)
            #print ("traverse weights: ", traverse_weights)


            # first assign a residence time to each node that is valid for this
            # path sampling iteration
            for node_keys in traverse_nodes:
                node = self.get_node(node_keys)
                _ = node.set_residence_time(
                                        radius_samp,
                                        self.isotope_decay_constant,
                                        self.tank_time_default,
                                        self.tank_cyl_time_default,
                                        self.pipe_time_default,
                                        )


            # calculate the inlet activity of the first node
            inlet_node = self.get_node(traverse_nodes[0])
            act_in = inlet_node.get_external_inflow_activity()
            res_time_in = 0.0
            for node_keys  in traverse_nodes:
                node = self.get_node(node_keys)
                act_in,res_time_in, = node.mc_node_calculation(
                                               act_in,
                                               res_time_in,
                                               self.tank_time_default,
                                               self.tank_cyl_time_default,
                                               self.pipe_time_default,
                                               self.isotope_decay_constant,
                                               res_t_fraction = 1,
                                               steady_state = True,
                                               outer_loop = True,
                                               )

            i += 1

            if (i) % step == 0:
                check = self.eval_mc_checkpoint(self.mc_error_max,
                                                self.n_sampling_min,i)
                if check:
                    break

            if i >= self.n_sampling_max:
                break

        self.write_cdgs_single_isotope()
        self.write_vtk_single_isotope()
        self.print_solution_single_isotope_steady_state()

        return 0

    def eval_mc_checkpoint(self,error, max_it,itn):
        """
        this function check in all points which is the number of iterations
        and the sampling error. When the conditions are met in all the nodes
        it returns true
        """

        check = True

        sum_counter = 0
        sum_error = 0


        n_nodes = len(self.all_nodes)

        max_err = 0
        max_it = 1e100

        for node_keys in self.all_nodes:

            node = self.get_node(node_keys)
            node.eval_mc_checkpoint()


            node_counter = node.mc_counter

            node_error = node.mc_error


            if node_counter < max_it:
                max_it = node_counter

            if node_error > max_err:
                max_err = node_error

            sum_counter += node_counter
            sum_error += node_error

            if node_counter < max_it or node_error > error:
                check = False

        av_counter = sum_counter/n_nodes
        av_error = sum_error/n_nodes


        it_string = f"iteration: {itn:>1.3e}"

        count_string = f"average number of iterations: {av_counter:>8.1f}"

        av_string =    f"average rel error: {av_error:>6.3f}"

        max_string =    f"max rel error: {max_err:>6.3f}"

        out_string = (it_string + "  "
                      + count_string + "  "
                      + av_string + "  "
                        + max_string)

        print (out_string)


        return check


    def finalize_mc_solve_single_isotope(self):
        """
        this function finalizes the monte carlo problem for a single isotope
        in steady state
        """
        for node_key in self.all_nodes:

            node = self.get_node(node_key)
            node.eval_mc_checkpoint()

        return 0


    def reset_node_counters(self):
        """
        this function reset the monte carlo counters for all the nodes
        """
        for node_keys in self.all_nodes:
            node = self.get_node(node_keys)
            node.reset_mc_counters(self.n_sources)

        return 0


    def print_solution_single_isotope_transient(self,time):
        """
        this function prints the tabular output of the transient simulation
        """
        header_cols = [
                        "circuit",
                        "group_name",
                        "node_id",
                        "length [m]",
                        "inlet_activity [Bq/m3]",
                        "average_activity [Bq/m3]",
                        "total_activity [Bq]",
                        "mc_relative_error",
                         ]

        header_string = "{: >30},"*len(header_cols) + "\n"

        data_string = ("{: >30},{: >30},{: >30},{: >30.6e},{: >30.6e},{: >30.6e},"
                       "{: >30.6e},{: >30.6e}\n")


        for circuit in self.circuits.values():

            out_name = f"_output_transient_{time:.3f}_sec.csv"
            file_name = circuit.name + out_name
            file_path = os.path.join(circuit.full_path,file_name)



            with open(file_path, 'w',encoding='utf-8') as f_out:
                f_out.write(header_string.format(*header_cols))


                for node in circuit.circuit_dictionary.values():
                    f_out.write(data_string.format(
                                                circuit.name,
                                                node.group_name,
                                                node.node_id,
                                                node.length_cm/100,
                                                node.inlet_activity_bq_m3,
                                                node.average_activity_bq_m3,
                                                node.tot_activity_bq,
                                                node.mc_error,
                                                ))




        return 0
    def print_solution_single_isotope_steady_state(self):
        """
        this function prints the tabular output of the steady state simulation
        """
        header_cols = [
                        "circuit",
                        "group_name",
                        "node_id",
                        "node_type",
                        "length [m]",
                        "node_volume [m3]",
                        "average_activity [Bq/m3]",
                        "outlet_activity [Bq/m3]",
                        "total_activity [Bq]",
                        "average_residence_time [s]",
                        "cumulated_residence_time [s]",
                        "linear_velocity [m/s]",
                        "density [kg/m3]",
                        "mass_flow [kg/s]",
                        "mc_relative_error",
                        "mc_sampling_n",
                         ]

        header_string = "{: >30},"*len(header_cols) + "\n"

        data_string=("{: >30},{: >30},{: >30},{: >30},{: >30.6e},{: >30.6e},{: >30.6e},{: >30.6e},"
                     "{: >30.6e},{: >30.6e},{: >30.6e},{: >30.6e},{: >30.6e},"
                     "{: >30.6e},{: >30.6e},{: >30d}\n")


        for circuit in self.circuits.values():

            out_name = "_output_steady_state.csv"
            file_name = circuit.name + out_name
            file_path = os.path.join(circuit.full_path,file_name)


            with open(file_path, 'w',encoding='utf-8') as f_out:
                f_out.write(header_string.format(*header_cols))


                for node in circuit.circuit_dictionary.values():
                    f_out.write(data_string.format(
                                                circuit.name,
                                                node.group_name,
                                                node.node_id,
                                                node.node_type,
                                                node.length_cm/100,
                                                node.volume_cm3/1e6,
                                                node.average_activity_bq_m3,
                                                node.mc_out_activity_bq_m3,
                                                node.tot_activity_bq,
                                                node.mc_res_time,
                                                node.mc_out_time_cumulated,
                                                node.linear_velocity_m_s,
                                                node.density_g_cm3*1000,
                                                node.mass_flow,
                                                node.mc_error,
                                                node.mc_counter,
                                                ))




        return 0




    def write_cdgs_single_isotope(self,t_step=''):
        """
        this function write the cdgs for a single isotope
        """


        for circuit in self.circuits.values():

            if t_step != '':
                file_name = circuit.name + "_" +  str(int(t_step)) +"sec.cdgs"
            else:
                file_name = circuit.name + ".cdgs"
            file_path = os.path.join(circuit.full_path,file_name)

            make_pipes_cdgs(circuit.circuit_dictionary,
                            self.isotope_decay_constant,
                            self.isotope_branching_ratio,
                            self.isotope_energy_bins,
                            self.isotope_energy_probs,
                            file_path,
                            )


        return 0


    def write_vtk_single_isotope(self,step=''):
        """
        this function write the vtk for a single isotope for transients
        the scale of the file is in meters
        """


        for circuit in self.circuits.values():

            if step != '':
                file_name = circuit.name + "_meters_tstep_" + step
            else:
                file_name = circuit.name + "_meters"

            file_path = os.path.join(circuit.full_path,file_name)

            write_vtk_files(circuit.circuit_dictionary,
                                file_path,
                                self.isotope_decay_constant,
                                steady_state_flag=self.steady_state,
                                )


        return 0


    def get_inflow_node_list(self):
        """
        this function returns a list of all the source nodes
        """

        source_node_list = []


        for node in self.all_nodes:
            node_test = self.get_node(node)
            inflow = node_test.get_external_inflow()
            if inflow != 0 :
                source_node_list.append(node)

        return source_node_list

    def get_outflow_node_list(self):
        """
        this function returns a list of all the loose nodes where fluid
        is leaving the circuit
        """

        loose_node_list = []


        for node in self.all_nodes:
            node_test = self.get_node(node)
            outflow = node_test.get_external_outflow()
            if outflow != 0:
                loose_node_list.append(node)


        return loose_node_list

    def get_total_outflow(self):
        """
        this function calculates the total outflow of the circuit
        """

        total_outflow = 0


        for node in self.all_nodes:
            sel_node = self.get_node(node)
            total_outflow += sel_node.get_external_outflow()

        return total_outflow

    def get_total_inflow(self):
        """
        this function returns the sum of all the inflows in the circuit
        """

        inflow_sum = 0

        for node in self.all_nodes:
            node_test = self.get_node(node)
            inflow_sum += node_test.get_external_inflow()

        return inflow_sum


    def traverse_nodes_forward(self, node):
        """
        this function takes the couple of keys that identify a node:
        circuit key + node key and, from there walks the chain until it
        reaches loose nodes
        """

        walked_nodes_list = np.empty(self.tot_n_nodes, dtype=object)
        walked_nodes_weights = np.ones(self.tot_n_nodes, dtype=object)

        walked_nodes_list[0] = node

        i = 1

        while True:

            current_node = self.get_node(node)

            child_node_keys,ds_weight = current_node.random_weighted_child()


            if not child_node_keys:
                break


            node = child_node_keys

            walked_nodes_list[i] = node
            walked_nodes_weights[i] = ds_weight
            i += 1

        walked_nodes_list = walked_nodes_list[:i]
        walked_nodes_weights = walked_nodes_weights[:i]


        return walked_nodes_list,walked_nodes_weights

    def traverse_nodes_forward_unweighted(self, node):
        """
        this function takes the couple of keys that identify a node:
        circuit key + node key and, from there walks the chain until it
        reaches loose nodes. the child node is chosen randomly
        """

        walked_nodes_list = np.empty(self.tot_n_nodes, dtype=object)
        walked_nodes_weights = np.ones(self.tot_n_nodes, dtype=object)

        walked_nodes_list[0] = node

        i = 1

        while True:

            current_node = self.get_node(node)

            child_node_keys,ds_weight = current_node.random_unweighted_child()


            if not child_node_keys:
                break


            node = child_node_keys

            walked_nodes_list[i] = node
            walked_nodes_weights[i] = ds_weight
            i += 1

        walked_nodes_list = walked_nodes_list[:i]
        walked_nodes_weights = walked_nodes_weights[:i]


        return walked_nodes_list,walked_nodes_weights

    def get_node(self, node_keys):
        """
        this function takes as input a tuple of node keys and
        returns the node object of the circuit
        """

        try:
            selected_node = (
                  self.circuits[node_keys[0]].circuit_dictionary[node_keys[1]]
                             )
        except KeyError:
            selected_node = False

        return selected_node

    def define_start_nodes(self):
        """
        this function returns the list of the possible start nodes and the
        corresponding probabilities. The list of the starting nodes contains also
        an identifier for the source. All the source-nodes, if present are
        grouped with the first identifier. All the nodes with a reaction rate
        have their own identifier. The only exception is for those nodes that
        are linked to an source-node and have also a reaction rate.
        """

        ext_source_node_list = []

        ext_source_node_probs = []

        irradiated_node_list = []

        irradiated_node_probs = []


        start_node_list = []

        start_node_probs = []

        counter_list = []



        mc_id_counter = 0
        irr_start_list = [0]


        for node in self.all_nodes:
            node_test = self.get_node(node)
            inflow = node_test.get_external_inflow()
            if inflow != 0 :
                ext_source_node_list.append(node_test.node_keys)
                ext_source_node_probs.append(node_test.mass_flow)

        for node in self.all_nodes:
            node_test = self.get_node(node)
            reac_rate = node_test.reaction_rate_m3

            #if reac_rate != 0 and node_test.node_keys not in ext_source_node_list:
            if reac_rate != 0:
                irradiated_node_list.append(node_test.node_keys)
                irradiated_node_probs.append(node_test.mass_flow)

        for node, mass_flow in zip(ext_source_node_list, ext_source_node_probs):
            start_dict = {}
            start_dict["node"] = node
            start_dict["mass_flow"] = mass_flow
            start_dict["counter_id"] = mc_id_counter
            sel_node = self.get_node(node)
            #sel_node.mc_id.append(mc_id_counter)
            start_node_list.append(start_dict)
            start_node_probs.append(mass_flow)
            #start_node_probs.append(1)
            #counter_list.append(mc_counter)
            #mc_counter += 1

        counter_list = [0]
        mc_id_counter += 1

        for node, mass_flow in zip(irradiated_node_list, irradiated_node_probs):
            #start_dict = {}
            #start_dict["node"] = node
            #start_dict["mass_flow"] = mass_flow
            #start_dict["counter_id"] = mc_counter

            irr_start_list.append(mc_id_counter)
            sel_node = self.get_node(node)
            sel_node.mc_id = mc_id_counter

            #start_node_list.append(start_dict)
            #start_node_probs.append(mass_flow)
            #start_node_probs.append(1)
            #start_node_list.append(start_dict)
            #start_node_probs.append(mass_flow)
            #start_node_probs.append(1)

            counter_list.append(mc_id_counter)

            mc_id_counter += 1



        tot_prob = sum(start_node_probs)
        start_node_probs[:] = [x / tot_prob for x in start_node_probs]


        #print ("start_node_list", start_node_list)
        #print ("start_node_probs", start_node_probs)
        #print("mc_counter_list", counter_list)
        #print ("irr_start_list", irr_start_list)

        return start_node_list, start_node_probs, irr_start_list, counter_list

    def pick_start_mc_counter(self):
        """
        this function returns a random choice of the monte carlo counters
        """
        mc_counter = np.random.choice(self.mc_list,1)[0]

        return mc_counter

    def pick_start_node(self):
        """
        this function returns a weighted random choice of the starting node
        INDEPENDENT COUNTERS - SINGLE SAMPLE POOL version
        """


        start_dict = np.random.choice(
                        self.ext_source_list,
                        1,
                        p=self.ext_source_prob,
                        )[0]
        start_node = start_dict['node']
        start_idx = start_dict['counter_id']
        if start_idx == 0:
            irr_flag = False
        else:
            irr_flag = True


        return start_node, start_idx, irr_flag




    def get_all_nodes(self):
        """
        this function return a vector of all the nodes of the case.
        the vector is a vector of tuple of node keys (circuit_name,node_id)
        """

        nodes_vector = []

        for circuit_key, circuit_def in self.circuits.items():
            for node_id in circuit_def.circuit_dictionary.keys():
                nodes_vector.append((circuit_key, node_id))

        return nodes_vector


    def get_all_attr(self, attr):
        """
        this function returns a dictionary with all the node values for
        a specific attribute
        """

        value_dict = {}


        for node in self.all_nodes:
            value_dict[node] = getattr(self.get_node(node), attr)

        return value_dict

    def print_all_nodes(self):
        """
        this function prints all the nodes of the case
        """
        for node_keys in self.all_nodes:
            node = self.get_node(node_keys)
            node.print_node()

        return 0
