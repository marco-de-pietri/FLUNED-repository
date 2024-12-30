"""
this file containes the class of the node object
"""
import math
import random
import sys
import numpy as np

from iapws.iapws97 import IAPWS97

from .utils import bin_search, linear_interp
from fluned.of_class import SimulationOF

from .single_isotope_activation import (
    average_activity,
    rr_average_activity,
    outlet_activity,
    rr_outlet_activity,
    res_time_from_reduction_rate,
)


class CircuitNode:
    """
    class to define a circuit node
    """

    def __init__(self, **node_params):
        """
        initialize a circuit node
        """

        self.node_id = node_params["node_id"]
        self.circuit_name = node_params["node_circuit"]
        self.node_keys = (self.circuit_name, self.node_id)
        self.parents = {}
        self.children = {}
        self.children_keys = []
        self.children_probs = []
        self.children_weights = []
        self.children_indexes = []
        self.node_type = node_params["node_type"]
        self.node_type_value = node_params["node_type_value"]
        self.activity_scaling = node_params["activity_scaling"]
        self.mass_flow = 0
        self.density_g_cm3 = 0
        self.kinematic_viscosity_m2_s = 0
        self.reynolds = 0
        self.n_turb = 0
        self.temperature_k = node_params["temperature_k"]
        self.pressure_mpa = node_params["pressure_mpa"]
        self.volume_cm3 = node_params["volume_cm3"]
        self.mcnp_cell = node_params["mcnp_cell"]
        self.axis_x = node_params["axis_x"]
        self.axis_y = node_params["axis_y"]
        self.axis_z = node_params["axis_z"]
        self.axis_x_1 = node_params["axis_x_1"]
        self.axis_y_1 = node_params["axis_y_1"]
        self.axis_z_1 = node_params["axis_z_1"]
        self.origin_x_cm = node_params["origin_x_cm"]
        self.origin_y_cm = node_params["origin_y_cm"]
        self.origin_z_cm = node_params["origin_z_cm"]
        self.length_cm = node_params["length_cm"]
        self.radius_cm = node_params["radius_cm"]
        self.group_name = node_params["group_name"]

        self.rtd_time = node_params["node_rtd_time"]
        self.rtd_data = node_params["node_rtd_data"]
        self.cfd_path = node_params["node_cfd_path"]

        if self.cfd_path != '':
            # initialize the cfd data
            self.cfd_sim = SimulationOF(self.cfd_path)
            self.cfd_sim.post_process_simulation()
            # update the volume with the data from the cfd simulation
            self.volume_cm3 = self.cfd_sim.volume_m3*1e6
        else:
            self.cfd_sim = None



        self.rtd_norm = []
        self.rtd_cumulative_activity_outlet = []
        self.rtd_cumulative_activity_average = []

        self.bulk_res_time = 0
        self.bulk_velocity_m_s = 0
        self.bulk_decay_rate_norm = 0
        self.bulk_average_activity_norm = 0
        self.u_max_norm = 0
        self.u_max_turb = 0
        self.u_max_laminar = 0


        self.sample_res_time = 0
        self.sample_a_frac = 0

        self.tot_activity_bq = 0
        self.reaction_rate_m3 = 1e6*node_params["reaction_rate_cm3"]

        self.mc_id = 0

        self.mc_error = 0

        self.mc_counter = 0
        self.mc_counter_vec = []



        self.mc_average_activity_bq_m3 = 0
        self.mc_average_activity_bq_m3_sum = 0
        self.mc_average_activity_bq_m3_vec = []

        self.mc_average_activity_bq_m3_sq = 0
        self.mc_average_activity_bq_m3_sq_sum = 0
        self.mc_average_activity_bq_m3_sq_vec = []

        self.inlet_activity_bq_m3 = 0

        self.mc_out_activity_bq_m3 = 0
        self.mc_out_activity_bq_m3_sum = 0
        self.mc_out_activity_bq_m3_vec = []

        self.mc_res_time = 0
        self.mc_res_time_sum = 0
        self.mc_res_time_vec = []

        self.mc_out_time_cumulated = 0
        self.mc_out_time_cumulated_sum = 0
        self.mc_out_time_cumulated_vec = []

        self.linear_velocity_m_s = 0

    def print_node(self):
        """
        print a node definition
        """

        print("node ID: ", self.node_id)
        print("node line name: ", self.group_name)
        print("node circuit name: ", self.circuit_name)
        print("node type           : ", self.node_type)
        print("mass Flow   [kg/s]  : ", self.mass_flow)
        print("Volume      [cm3]   : ", self.volume_cm3)
        print("Temperature [K]     : ", self.temperature_k)
        print("Pressure    [MPa]   : ", self.pressure_mpa)
        print("density     [g/cm3] : ", self.density_g_cm3)
        print("number of parents   : ", len(self.parents))
        for parent, parent_dict in self.parents.items():
            print("Parent ID: ", parent)
            print("\tParent type: ", parent_dict["parent_type"])
            print("\tParent external type: ", parent_dict["parent_ext_type"])
            print("\tParent external mass inflow: ", parent_dict["parent_ext_mass_inflow"])
            print("\tParent ext-node mass inflow: ", parent_dict["parent_ext_node_mass_inflow"])
            print("\tParent internal mass inflow: ", parent_dict["parent_int_mass_inflow"])
            print("\tParent node id: ", parent_dict["parent_node_id"])
            print("\tParent circuit: ", parent_dict["parent_circuit"])
            print("\tParent fraction: ", parent_dict["fraction"])
            print("\tParent keys: ", parent_dict["parent_keys"])

        print("number of children   : ", len(self.children))
        for child, child_dict in self.children.items():
            print("Child ID: ", child)
            print("\tChild type: ", child_dict["child_type"])
            print("\tChild external type: ", child_dict["child_ext_type"])
            print("\tChild node id: ", child_dict["child_node_id"])
            print("\tChild circuit: ", child_dict["child_circuit"])
            print("\tChild fraction: ", child_dict["fraction"])
            print("\tChild downstream fraction: ", child_dict["fraction_dstream"])
            print("\tChild keys: ", child_dict["child_keys"])
        print("\n")
        return 0

    def add_parent(self,
                   parent_vec,
                   parent_circuit = '',
                   ext_type = '',
                   ext_mass_inflow = 0,
                   ext_inflow_activity = 0,
                   ext_inflow_time = 0,
                   int_mass_inflow = 0,
                   ext_node_mass_inflow = 0,
                   ):
        """
        add a parent to the node definition
        """
        parent_dic = {}
        parent_dic["parent_id"] = int(parent_vec[0])
        parent_dic["parent_circuit"] = parent_circuit
        parent_dic["parent_type"] = parent_vec[1]
        parent_dic["parent_ext_type"] = ext_type
        parent_dic["parent_ext_mass_inflow"] = ext_mass_inflow
        parent_dic["parent_int_mass_inflow"] = int_mass_inflow
        parent_dic["parent_ext_node_mass_inflow"] = ext_node_mass_inflow
        parent_dic["parent_ext_activity_inflow"] = ext_inflow_activity
        parent_dic["parent_ext_activity_time"] = ext_inflow_time
        parent_dic["parent_node_id"] = int(parent_vec[2])
        parent_dic["fraction"] = float(parent_vec[3])
        parent_dic["parent_keys"]=(parent_circuit, parent_dic["parent_node_id"])

        self.parents[int(parent_vec[0])] = parent_dic

    def add_child(
        self,
        child_node_id,
        child_type,
        child_fraction,
        child_circuit,
        ext_type = '',
        child_fraction_dstream = 0,
    ):
        """
        add a children to the node definition
        """
        new_id = len(self.children) + 1
        child_dict = {}
        child_dict = {
            "child_node_id": child_node_id,
            "child_type": child_type,
            "child_ext_type": ext_type,
            "child_circuit": child_circuit,
            "fraction": child_fraction,
            "fraction_dstream": child_fraction_dstream,
            "child_keys": (child_circuit, child_node_id),
        }

        self.children[new_id] = child_dict

        return 0

    def get_external_inflow(self):
        """
        this function returns the external inflow of the node coming from
        sources of source-node type
        """

        ext_inflow_mass = 0


        for parent in self.parents.values():
            ext_inflow_mass += parent["parent_ext_mass_inflow"]


        return ext_inflow_mass

    def get_external_inflow_activity(self,t_sample = 0, mc_counter = 0):
        """
        This function returns the amount activity entering into the node from
        external source (source-node type)
        """

        if  mc_counter != 0:
            return 0


        temp_mass = 0
        temp_act = 0

        for parent in self.parents.values():

            if ((parent["parent_type"] == "ext" and
                parent["parent_ext_type"] == "source")):

                if parent["parent_ext_activity_inflow"] != -2:
                    temp_mass += parent["parent_ext_mass_inflow"]
                    temp_act+= (parent["parent_ext_mass_inflow"]*
                                parent["parent_ext_activity_inflow"])
                elif parent["parent_ext_activity_inflow"] == -2:
                    temp_mass += parent["parent_ext_mass_inflow"]
                    closest_left_idx = bin_search(
                                parent["parent_ext_activity_time"],
                                t_sample)
                    act_sample = linear_interp(
                        parent["parent_ext_activity_time"],
                        parent["parent_ext_activity_data"],
                        closest_left_idx,
                        t_sample,
                        outside = 'zero',
                     )

                    temp_act+= parent["parent_ext_mass_inflow"]*act_sample


        if temp_mass == 0:
            # this option is for nodes with activation
            act = 0
        else:
            act = temp_act/temp_mass

        return act


    def get_external_outflow(self):
        """
        this function returns the flow exiting the system
        """

        ext_outflow_mass = 0

        if len(self.children) == 0:
            ext_outflow_mass = self.mass_flow


        return ext_outflow_mass

    def define_children_weights(self):
        """
        this function defines the weights probabilities of the water to go
        into the different children nodes
        """

        child_keys = []
        child_probs = []
        child_weights = []
        child_indexes = []
        for i,child in enumerate(self.children.values()):
            child_keys.append(child["child_keys"])
            child_probs.append(child["fraction"])
            child_weights.append(child["fraction_dstream"])
            child_indexes.append(i)
        self.children_keys = child_keys
        self.children_probs = child_probs
        self.children_weights = child_weights
        self.children_indexes = child_indexes

        return 0

    def assign_node_water_properties(self):
        """
        this function assign a density to the node assuming it is water
        using the temperature and pressure. it uses the iapws module which
        implements the iapws-97 formulation
        """

        w_properties = IAPWS97(T=self.temperature_k, P=self.pressure_mpa)

        if w_properties.region != 1:
            print("water not in liquid region")
            print("temperature: [K]", self.temperature_k)
            print("pressure: [MPa]", self.pressure_mpa)
            sys.exit()

        self.density_g_cm3 = w_properties.rho / 1000

        self.kinematic_viscosity_m2_s = w_properties.nu

        return 0

    def calculate_bulk_properties(self):
        """
        once a density and a mass flow is assigned to each node we can
        calculate its bulk properties:

        reynolds
        n_turb
        bulk velocity
        uniform velocity - residence time

        u_max_norm
        u_max_turb
        u_max_laminar

        """

        node_vol_m3 = self.volume_cm3 * 1e-6
        density_kg_m3 = self.density_g_cm3 * 1e3
        vol_flow_m3_s = self.mass_flow / density_kg_m3
        length_m = self.length_cm * 1e-2
        diameter_m = self.radius_cm * 2 * 1e-2

        self.bulk_res_time = node_vol_m3 / vol_flow_m3_s

        if length_m != 0:
            self.bulk_velocity_m_s = length_m / self.bulk_res_time

            self.reynolds = (
                self.bulk_velocity_m_s * diameter_m /
                self.kinematic_viscosity_m2_s
                            )
            self.n_turb = (2 / 3) * math.log(self.reynolds)

            self.u_max_norm = (
                (self.n_turb + 1) *
                (2 * self.n_turb + 1) / (2 * (self.n_turb**2))
                              )

            self.u_max_turb = self.u_max_norm * self.bulk_velocity_m_s
            self.u_max_laminar = 2 * self.bulk_velocity_m_s

        return 0

    def calculate_bulk_activation(self, dec_constant):
        """
        this function calculates the bulk activation - decay of a node using
        the uniform velocity approximation. At the moment is limited to
        decay only
        uniform velocity - decay rate
        uniform velocity - normalized average activity
        """

        self.bulk_decay_rate_norm = outlet_activity(
            a_in_vol=1,
            rr_vol=0,
            lambda_decay=dec_constant,
            res_time=self.bulk_res_time,
        )

        self.bulk_average_activity_norm = average_activity(
            a_in_vol=1,
            rr_vol=0,
            lambda_decay=dec_constant,
            res_time=self.bulk_res_time,
        )

        return 0

    def preprocess_rtd_data(self, dec_constant):
        """
        this function preprocess the rtd and creates arrays that can be used
        for the decay in components for which a residence time distribution is
        defined
        """
        if len(self.rtd_data) == 0:
            return 0

        norm_factor = sum(self.rtd_data)


        self.rtd_norm = self.rtd_data/norm_factor

        decay_array = np.zeros(len(self.rtd_norm))
        average_array = np.zeros(len(self.rtd_norm))


        for i, (t_frac, frac) in enumerate(zip(self.rtd_time, self.rtd_norm)):
            if i == 0:
                decay_array[i] = frac
                average_array[i] = frac
                continue
            decay_array[i] = frac*outlet_activity(
                                    a_in_vol=1,
                                    rr_vol=0,
                                    lambda_decay=dec_constant,
                                    res_time=t_frac,
                                     )
            average_array[i] = frac*average_activity(
                                   a_in_vol=1,
                                   rr_vol=0,
                                   lambda_decay=dec_constant,
                                   res_time=t_frac,
                                    )


        self.rtd_cumulative_activity_outlet = np.cumsum(decay_array)
        self.rtd_cumulative_activity_average = np.cumsum(average_array)

        return 0


    def random_unweighted_parent(self):
        """
        this function returns a casual parent of the node
        """

        random_parent = random.choice(list(self.parents.values()))

        return random_parent

    def random_weighted_child(self):
        """
        this function returns sample a child according to the respective
        properties
        """
        if len(self.children) == 1:
            for child in self.children.values():
                return child["child_keys"],child['fraction_dstream']

        try:
            child_idx = random.choices(
                self.children_indexes, weights=self.children_probs, k=1
            )[0]

            return self.children_keys[child_idx], self.children_weights[child_idx]

        #try:
        #    child_idx = random.choices(
        #        self.children_indexes,  k=1
        #    )[0]

        #    return self.children_keys[child_idx], self.children_weights[child_idx]

        except IndexError:
            return None, None

    def random_unweighted_child(self):
        """
        this function returns sample a child according to the respective
        properties
        """
        if len(self.children) == 1:
            for child in self.children.values():
                return child["child_keys"],child['fraction_dstream']


        try:
            child_idx = random.choices(
                self.children_indexes,  k=1
            )[0]

            return self.children_keys[child_idx], self.children_weights[child_idx]

        except IndexError:

            return None, None

    def node_calculation(
        self,
        act_in,
        t_in,
        tank_method,
        tank_cyl_method,
        pipe_method,
        lambda_decay,
        res_t_fraction=1,
        res_t_complete_fraction = 1,
        steady_state = True,
        outer_loop = True,
        mc_counter = 0,
    ):
        """
        this function manages the calculation of the activity and the
        residence time of fluid crossing a node.

        first, it calls the function to calculate the average and outlet activity,
        get_outlet_average_activity

        Finally it updates the monte carlo tallies and returns the outlet time
        and activity.

        """

        lower_limit = 1e-20



        # this flag allows to activate the node activation only in one
        # source at the time


        if steady_state:
            activ_flag = True
            mc_counter = 0

        if not steady_state:
            if mc_counter == self.mc_id and mc_counter != 0:
                activ_flag = True
            else:
                activ_flag = False

        # outer loops
        if outer_loop:
            raise_tallies = True

        # inner loops
        if not outer_loop:
            raise_tallies = False






        act_out, act_average = self.get_outlet_average_activity(
            act_in,
            tank_method,
            tank_cyl_method,
            pipe_method,
            lambda_decay,
            res_t_fraction = res_t_fraction,
            res_t_complete_fraction = res_t_complete_fraction,
            steady_state = steady_state,
            activation_flag = activ_flag,
        )


        # check if the activity is below the lower limit
        if act_out < lower_limit:
            act_out = 0


        res_time = self.sample_res_time * res_t_fraction
        t_out = t_in + res_time

        # update the monte carlo tallies
        if raise_tallies:

            #if self.node_id == 10:


            #    #print ("act_in initial", act_in_1)
            #    #print ("weight", weight)
            #    print ("mc_counter", mc_counter)
            #    print ("act_in", act_in)
            #    print ("act_out", act_out)

            #self.mc_average_activity_bq_m3_sum += act_average*weight
            #self.mc_average_activity_bq_m3_vec[mc_counter] += act_average*weight

            #self.mc_average_activity_bq_m3_sq_sum += ((act_average*weight)**2)
            #self.mc_average_activity_bq_m3_sq_vec[mc_counter] += ((act_average*weight)**2)

            #self.mc_out_activity_bq_m3_sum += act_out*weight
            #self.mc_out_activity_bq_m3_vec[mc_counter] += act_out*weight

            #self.mc_res_time_sum += self.sample_res_time * res_t_fraction*weight
            #self.mc_res_time_vec[mc_counter] += self.sample_res_time * res_t_fraction*weight

            #self.mc_out_time_cumulated_sum += t_out*weight
            #self.mc_out_time_cumulated_vec[mc_counter] += t_out*weight

            self.mc_average_activity_bq_m3_sum += act_average
            self.mc_average_activity_bq_m3_vec[mc_counter] += act_average

            self.mc_average_activity_bq_m3_sq_sum += ((act_average)**2)
            self.mc_average_activity_bq_m3_sq_vec[mc_counter] += ((act_average)**2)

            self.mc_out_activity_bq_m3_sum += act_out
            self.mc_out_activity_bq_m3_vec[mc_counter] += act_out

            self.mc_res_time_sum += res_time
            self.mc_res_time_vec[mc_counter] += res_time

            self.mc_out_time_cumulated_sum += t_out
            self.mc_out_time_cumulated_vec[mc_counter] += t_out


            self.increase_mc_counter(mc_counter)

        return act_out, t_out

    def get_outlet_average_activity(
        self,
        act_in,
        tank_method,
        tank_cyl_method,
        pipe_method,
        lambda_decay,
        res_t_fraction = 1,
        res_t_complete_fraction = 1,
        steady_state = True,
        activation_flag = False,
    ):
        """
        this function calculates the average and outlet activity when a fluid
        crosses the node. It takes as input the residence time calculated
        with set_residence_time.

        if the method value is the default (-1), it calls the default method
        for the node type

        if the method value corresponds to a reduction rate (between -1 and 0),
        it uses a back-calculated residence time from the reduction rate. if the
        reaction is different than zero it gives an error.

        if the method value corresponds to a user defined residence time (
        greater than 0), the outlet activity is calculated accordingly
        """

        act_out = 0
        act_average = 0


        if activation_flag:
            reac_rate = self.reaction_rate_m3
        else:
            reac_rate = 0

        # default methods
        if self.node_type_value == -1:
            if self.node_type == "tank":
                method = tank_method
                act_out, act_average = self.calculate_outlet_average_activity_method(
                    act_in,
                    lambda_decay,
                    self.sample_res_time,
                    reac_rate,
                    method,
                    res_t_fraction,
                    res_t_complete_fraction,
                    steady_state,
                )

            elif self.node_type == "tank-cyl":
                method = tank_cyl_method
                act_out, act_average = self.calculate_outlet_average_activity_method(
                    act_in,
                    lambda_decay,
                    self.sample_res_time,
                    reac_rate,
                    method,
                    res_t_fraction,
                    res_t_complete_fraction,
                    steady_state,
                )

            elif self.node_type == "pipe":
                method = pipe_method
                act_out, act_average = self.calculate_outlet_average_activity_method(
                    act_in,
                    lambda_decay,
                    self.sample_res_time,
                    reac_rate,
                    method,
                    res_t_fraction,
                    res_t_complete_fraction,
                    steady_state,
                )

            else:
                print("ERROR: node type not recognized")
                print(self.node_type)
                sys.exit()

        # reduction rate method (not for activation)
        elif -1 < self.node_type_value < 0:
            if self.reaction_rate_m3 != 0:
                print("ERROR: red-rate method not applicable for activation")
                print("ERROR NODE: ",self.node_keys)
                sys.exit()

            res_time_calc = self.sample_res_time*res_t_fraction

            act_average = average_activity(
                a_in_vol=act_in,
                rr_vol=0,
                lambda_decay=lambda_decay,
                res_time=res_time_calc,
            )
            act_out = -self.node_type_value * act_in

            if res_t_fraction != 1:
                act_out = 0
                act_average = act_average * res_t_fraction

        # no decay method
        elif self.node_type_value == 0:

            if self.reaction_rate_m3 != 0:
                print("ERROR: no-decay method not applicable for activation")
                print("ERROR NODE: ",self.node_keys)
                print(self.node_type_value)
                sys.exit()

            act_out = act_in
            act_average = act_in

        # user defined residence time method
        elif self.node_type_value > 0:

            if steady_state:
                res_time_calc = self.sample_res_time
                act_out = outlet_activity(
                    a_in_vol=act_in,
                    rr_vol=reac_rate,
                    lambda_decay=lambda_decay,
                    res_time=res_time_calc,
                )

                act_average = average_activity(
                    a_in_vol=act_in,
                    rr_vol=reac_rate,
                    lambda_decay=lambda_decay,
                    res_time = res_time_calc,
                )

            if not steady_state:

                #copy the implementation from the uniform method
                print ("ERROR: not steady state not implemented yet")
                sys.exit()

                # contribution from decay
                #res_time_decay = self.sample_res_time*res_t_fraction

                #if res_t_fraction < 1:
                #    dec_out = 0
                #else:
                #    dec_out = outlet_activity(
                #        a_in_vol=act_in,
                #        rr_vol=0,
                #        lambda_decay=lambda_decay,
                #        res_time=res_time_decay,
                #        )

                #dec_average = average_activity(
                #    a_in_vol=act_in,
                #    rr_vol=0,
                #    lambda_decay=lambda_decay,
                #    res_time = res_time_decay,
                #)

                #if reac_rate != 0:
                #    res_time_activ = self.sample_res_time*res_a_fraction

                #    if res_time_activ > res_time_decay:
                #        res_time_activ = np.random.uniform(0,res_time_decay)
                #        res_a_fraction = res_time_activ/self.sample_res_time

                #    rr_out = outlet_activity(
                #        a_in_vol=0,
                #        rr_vol=reac_rate,
                #        lambda_decay=lambda_decay,
                #        res_time=res_time_activ,
                #        )

                #    rr_average = average_activity(
                #         a_in_vol=0,
                #         rr_vol=reac_rate,
                #         lambda_decay=lambda_decay,
                #         res_time=res_time_activ,
                #         )

                #else:
                #    rr_out = 0
                #    rr_average = 0

                #act_out = rr_out + dec_out
                #act_average = ((dec_average*res_t_fraction) +
                #              (rr_average*res_a_fraction))

        elif self.node_type_value == -2:
            if self.node_type in ["tank", "tank-cyl","pipe"]:
                if steady_state:

                    dec_out = act_in * self.rtd_cumulative_activity_outlet[-1]
                    dec_average = act_in * self.rtd_cumulative_activity_average[-1]

                    if reac_rate == 0:

                        rr_out = 0
                        rr_average = 0

                    elif reac_rate > 0:


                        res_time_calc = self.sample_res_time
                        rr_out = outlet_activity(
                            a_in_vol=act_in,
                            rr_vol=reac_rate,
                            lambda_decay=lambda_decay,
                            res_time=res_time_calc,
                        )

                        rr_average = average_activity(
                            a_in_vol=act_in,
                            rr_vol=reac_rate,
                            lambda_decay=lambda_decay,
                            res_time=res_time_calc,
                        )

                    act_out = rr_out + dec_out
                    act_average = dec_average + rr_average

                # not steady state
                elif not steady_state:


                    # we are considering that the parcel sampled through the rtd
                    # goes through the node, no matter how fast

                    # res_t_fraction is used to consider if the residence time we
                    # sampled is longer than the residual time to reach the
                    # timestamp we are simulating (t_stamp - t_in < res_time)




                    res_time_decay = self.sample_res_time*res_t_fraction
                    res_a_fraction = self.sample_a_frac

                    if res_t_fraction < 1 or res_a_fraction < 1:
                        dec_out = 0

                    else:
                        dec_out = outlet_activity(
                            a_in_vol=act_in,
                            rr_vol=reac_rate,
                            lambda_decay=lambda_decay,
                            res_time=res_time_decay,
                        )

                    dec_average = average_activity(
                        a_in_vol=act_in,
                        rr_vol=reac_rate,
                        lambda_decay=lambda_decay,
                        res_time = res_time_decay,
                    )

                    # no activation
                    if reac_rate == 0:

                        rr_out = 0
                        rr_average = 0

                    # constant activation value
                    elif reac_rate > 0:

                        res_time_activ = self.sample_res_time*res_t_complete_fraction

                        rr_integrate_factor = self.sample_res_time*res_t_complete_fraction/res_a_fraction


                        rr_out = rr_outlet_activity(
                            rr_vol=reac_rate,
                            lambda_decay=lambda_decay,
                            res_time=res_time_activ,
                            )

                        rr_out = rr_out*rr_integrate_factor


                        rr_average = rr_average_activity(
                            rr_vol=reac_rate,
                            lambda_decay=lambda_decay,
                            res_time=res_time_activ,
                            )
                        rr_average = rr_average*rr_integrate_factor




                    #if self.node_id == 1 and rr_out != 0:
                    #    print ("rr_out", rr_out)
                    #    print ("dec_out", dec_out)
                    #    print ("res time", res_time)
                    #    print ("res time activ", res_time_activ)


                    act_out = rr_out + dec_out
                    act_average = ((dec_average*res_t_fraction) +
                                  (rr_average*res_t_complete_fraction))


                    #    #print("new time and interp")
                    #    #print (self.sample_res_time)
                    #    print(np.interp(self.sample_res_time,
                    #                        self.rtd_time,
                    #                        self.rtd_cumulative_activity_outlet,
                    #                        ))

                    #    act_out = act_in*np.interp(self.sample_res_time,
                    #                        self.rtd_time,
                    #                        self.rtd_cumulative_activity_outlet,
                    #                        )

                    ##act_average = act_in * self.rtd_cumulative_activity_average[idx_res_time]

            if self.node_type in ["cfd"]:
                if  steady_state:

                    act_out = act_in * self.cfd_sim.reduction_rate[-1]
                    act_average = act_in * self.cfd_sim.normalized_average_td[-1]

                    #print ("act_out", act_out)
                    #print ("act_average", act_average)

                    if activation_flag:

                        act_out += lambda_decay*self.cfd_sim.outlet_rr_conc_atoms_m3[-1]
                        act_average += lambda_decay*self.cfd_sim.average_ta[-1]
                        #print ("act_out with rr", act_out)
                        #print ("act_average with rr", act_average)



                elif not steady_state:

                    #copy the implementation from the uniform method
                    raise NotImplementedError("cfd transient not implemented yet")


        else:
            print("ERROR: node type value not recognized")
            print(self.node_type_value)
            sys.exit()

        return act_out, act_average

    def calculate_outlet_average_activity_method(
        self,
        act_in,
        lambda_decay,
        res_time,
        reac_rate,
        method,
        res_t_fraction = 1,
        res_t_complete_fraction =1,
        steady_state = False,
    ):
        """
        this function calculates the activation in pipe and tank nodes using the
        appropriate method. the iteration sampe rate has already been computed
        and passed as input
        """


        if method == "reynolds":

            # need to use the sampled residence time according to the reynold
            # number
            if steady_state:
                act_out = outlet_activity(
                    a_in_vol=act_in,
                    rr_vol=reac_rate,
                    lambda_decay=lambda_decay,
                    res_time=res_time,
                )

                act_average = average_activity(
                    a_in_vol=act_in,
                    rr_vol=reac_rate,
                    lambda_decay=lambda_decay,
                    res_time=res_time,
                )
            if not steady_state:

                #copy the implementation from the uniform method
                print ("ERROR not implemented yet")
                sys.exit()


        elif method == "uniform":
            if steady_state:
                if reac_rate == 0:
                    act_out = act_in * self.bulk_decay_rate_norm
                    act_average = act_in * self.bulk_average_activity_norm
                elif reac_rate != 0:
                    act_out = outlet_activity(
                        a_in_vol=act_in,
                        rr_vol=reac_rate,
                        lambda_decay=lambda_decay,
                        res_time=res_time,
                    )

                    act_average = average_activity(
                        a_in_vol=act_in,
                        rr_vol=reac_rate,
                        lambda_decay=lambda_decay,
                        res_time=res_time,
                    )
            if not steady_state:

                # res_t_fraction is used to consider if the residence time we
                # sampled is longer than the residual time to reach the
                # timestamp we are simulating (t_stamp - t_in < res_time)
                res_time_decay = res_time*res_t_fraction

                res_a_fraction = self.sample_a_frac

                #if self.node_id ==  1:
                #    print ("node id", self.node_id)
                #    print ("mc id", self.mc_id)
                #    print ("res_time", res_time)
                #    print ("res_t_fraction", res_t_fraction)
                #    print ("res_a_fraction", res_a_fraction)
                #    print ("react rate", reac_rate)


                # track does not reach the node outlet

                # if res_t_fraction < 1 means that the parcel does not reach the
                # node outlet

                # if res_a_fraction < 1 means that the we are tracking the reaction
                # rate ??

                if res_t_fraction < 1 or res_a_fraction < 1:
                    dec_out = 0

                # track reaches the node outlet
                else:
                    dec_out = outlet_activity(
                        a_in_vol=act_in,
                        rr_vol=0,
                        lambda_decay=lambda_decay,
                        res_time=res_time_decay,
                        )

                # the contribution to the average is computed regardless of
                # whether the track reaches the outlet or not
                dec_average = average_activity(
                    a_in_vol=act_in,
                    rr_vol=0,
                    lambda_decay=lambda_decay,
                    res_time = res_time_decay,
                )

                # outlet and average contribution due to the reaction rate
                if reac_rate != 0:


                    res_time_activ = res_time*res_t_complete_fraction

                    rr_integrate_factor = res_time*res_t_complete_fraction/res_a_fraction


                    #rr_t_contrib = reac_rate*lambda_decay*rr_integrate_factor

                    #rr_out = outlet_activity( a_in_vol=rr_t_contrib,
                    #    rr_vol=0,
                    #    lambda_decay=lambda_decay,
                    #    res_time=res_time_activ,
                    #    )

                    rr_out = rr_outlet_activity(
                        rr_vol=reac_rate,
                        lambda_decay=lambda_decay,
                        res_time=res_time_activ,
                        )

                    rr_out = rr_out*rr_integrate_factor

                    #rr_average = ((reac_rate/(lambda_decay*(res_time_activ**2))) *
                    #                (1-np.exp(-lambda_decay*res_time_activ)*
                    #                (1+lambda_decay*res_time_activ)))

                    rr_average = rr_average_activity(
                        rr_vol=reac_rate,
                        lambda_decay=lambda_decay,
                        res_time=res_time_activ,
                        )
                    rr_average = rr_average*rr_integrate_factor

                    if self.node_id == 1:
                        print ("node id", self.node_id)
                        print ("res_time", res_time)
                        print ("res_t_fraction", res_t_fraction)
                        print ("res_t_complete_fraction", res_t_complete_fraction)
                        print ("res_a_fraction", res_a_fraction)
                        print ("rr_integrate_factor", rr_integrate_factor)
                        print ("res_time_activ", res_time_activ)


                # no activation
                else:
                    rr_out = 0
                    rr_average = 0


                #if self.node_id == 1 and rr_out != 0:
                #    print ("rr_out", rr_out)
                #    print ("dec_out", dec_out)
                #    print ("res time", res_time)
                #    print ("res time activ", res_time_activ)


                act_out = rr_out + dec_out
                act_average = ((dec_average*res_t_fraction) +
                              (rr_average*res_t_complete_fraction))

        elif method == "no-decay":

            if self.reaction_rate_m3 != 0:
                print ("ERROR: no-decay method applied to node: ",self.node_keys)
                print ("with a reaction rate different than zero")
                sys.exit()

            act_out = act_in
            act_average = act_in

        else:
            print("ERROR: method not recognized")
            print(method)
            sys.exit()

        return act_out, act_average

    def set_residence_time(
        self,
        radius_norm,
        lambda_decay,
        tank_method,
        tank_cyl_method,
        pipe_method,
        found_source_flag = False,
        mc_counter = 0,
        numerical_method = 'monte_carlo'
    ):
        """
        this function calculates the time it takes to traverse a node.

        if the method value is for the value is the default (-1), it calls
        calculate_residence_time_methods which applies the deault methods

        if the method value corresponds to a reduction rate (between -1 and 0),
        it back-calculates the residence time from the reduction rate. if the
        reaction is different than zero it gives an error.

        if the method value corresponds to a user defined residence time (
        greater than 0), it uses the user defined residence time.
        """

        res_time = 0
        res_a_fraction = 1
        transient_activ_flag = False

        # this is used for the transient study of circuits with multiple
        # sources.
        if mc_counter == 0:
            found_source_flag = True
            transient_activ_flag = False
        else:
            if not found_source_flag:
                if mc_counter == self.mc_id:
                    transient_activ_flag = True
                    found_source_flag = True
                else:
                    self.sample_a_frac = res_a_fraction
                    self.sample_res_time = res_time
                    found_source_flag = False
                    return found_source_flag

        # default method
        if self.node_type_value == -1:
            if self.node_type == "tank":
                method = tank_method
                res_time = self.calculate_residence_time_method(
                    radius_norm,
                    method,
                )
            elif self.node_type == "tank-cyl":
                method = tank_cyl_method
                res_time = self.calculate_residence_time_method(
                    radius_norm,
                    method,
                )

            elif self.node_type == "pipe":
                method = pipe_method
                res_time = self.calculate_residence_time_method(
                    radius_norm,
                    method,
                )

            else:
                print("ERROR: node type not recognized")
                print(self.node_type)
                sys.exit()

            # if there is activation the residence time is reduced because
            # we assume that a particle can be created at any time during the
            # sampled path.
            if transient_activ_flag:
                res_a_fraction = np.random.uniform()
                res_time = res_time*res_a_fraction

        # reduction rate method (not for activation)
        elif -1 < self.node_type_value < 0:
            if self.reaction_rate_m3 != 0:
                print("ERROR: red-rate method not applicable for activation")
                print(self.node_type_value)
                sys.exit()

            res_time = res_time_from_reduction_rate(self.node_type_value,
                                                    lambda_decay)

        # no decay method
        elif self.node_type_value == 0:
            if self.reaction_rate_m3 != 0:
                print("ERROR: no-decay method not applicable for activation")
                print(self.node_type_value)
                sys.exit()

            res_time = 0


        # user defined residence time method
        elif self.node_type_value > 0:
            user_res_time = self.node_type_value
            res_time = user_res_time

            if transient_activ_flag:
                print ("ERROR not implemented yet")
                sys.exit()

        elif self.node_type_value == -2:
            # a rtd is defined
            if self.node_type in ["tank", "tank-cyl","pipe"]:
                res_time = random.choices(self.rtd_time, weights=self.rtd_norm)[0]

                if transient_activ_flag:
                    res_a_fraction = np.random.uniform()
                    res_time = res_time*res_a_fraction
            # FLUNED cfd simulation is provided
            elif self.node_type in "cfd":
                reduction_rate = self.cfd_sim.reduction_rate[-1]
                res_time = res_time_from_reduction_rate(reduction_rate,
                                                        lambda_decay)

        else:
            print("ERROR: node type value not recognized")
            print(self.node_type_value)
            raise ValueError


        self.sample_a_frac = res_a_fraction
        self.sample_res_time = res_time

        return found_source_flag

    def calculate_residence_time_method(
        self,
        radius_norm,
        method,
    ):
        """
        this function calculates the residence time in pipe and tank nodes
        using the appropriate method
        """

        if method == "reynolds":
            res_time = self.calculate_residence_time_reynolds(radius_norm)

        elif method == "uniform":
            res_time = self.bulk_res_time

        elif method == "no-decay":
            res_time = 0
        else:
            print("ERROR: method not recognized")
            print(method)
            raise ValueError

        return res_time

    def calculate_residence_time_reynolds(self, radius_norm):
        """
        this function calculates the residence time of a particle crossing
        a pipe along the adimensional radius radius_norm, depending using the
        raynolds number to calculate the a velocity profile
        """

        if self.reynolds <= 2300:
            vel = self.u_max_laminar * (1 - radius_norm**2)
            length_m = self.length_cm * 1e-2
            res_time = length_m / vel

        elif 2300 < self.reynolds < 4000:
            res_time = self.bulk_res_time

        elif self.reynolds >= 4000:
            vel = self.u_max_turb * (1 - radius_norm) ** (1 / self.n_turb)
            length_m = self.length_cm * 1e-2
            res_time = length_m / vel

        else:
            print("ERROR: Reynolds number not recognized")
            print(self.reynolds)
            sys.exit()

        return res_time

    def eval_mc_checkpoint(self):
        """
        this function calculates the average quantities for the node
        """

        self.calculate_mc_error_vec()
        self.calculate_mc_averages()


        return


    #def calculate_mc_error(self):
    #    """
    #    this function calculates the error of the monte carlo simulation
    #    the naive algorithm is used:
    #    https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    #    """

    #    mean = 0

    #    if self.mc_counter > 1:

    #        count = self.mc_counter
    #        act_sum = self.mc_average_activity_bq_m3_sum
    #        act_sq_sum = self.mc_average_activity_bq_m3_sq_sum


    #        mean = self.mc_average_activity_bq_m3_sum / self.mc_counter
    #        variance = (abs(act_sq_sum - (act_sum**2 / count))) / (count - 1)
    #        std_deviation_mean = math.sqrt(variance)/math.sqrt(count)


    #    if mean == 0:
    #        mc_error = 0
    #    else:
    #        mc_error = std_deviation_mean/mean

    #    self.mc_error = mc_error

    #    return 0

    def calculate_mc_error_vec(self):
        """
        this function calculates the error of the monte carlo simulation
        the naive algorithm is used:
        https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
        """

        node_abs_err = 0
        node_val = 0

        for i, _ in enumerate(self.mc_counter_vec):

            count = self.mc_counter_vec[i]

            if count > 1:
                act_sum = self.mc_average_activity_bq_m3_vec[i]
                act_sq_sum = self.mc_average_activity_bq_m3_sq_vec[i]
                mean = act_sum / count
                variance = (abs(act_sq_sum - (act_sum**2 / count))) / (count - 1)
                std_deviation_mean = math.sqrt(variance)/math.sqrt(count)

                node_abs_err += (std_deviation_mean**2)
                node_val += mean

        if node_val == 0:

            mc_error = 0

        else:

            mc_error = math.sqrt(node_abs_err)/node_val




        self.mc_error = mc_error

        return 0

    def write_cfd_vtk(self, file_path, decay_const):
        """
        this function writes the vtk file for the cfd node at the end of the
        calculation - it provides the inlet activity to scale the concentration
        field due only to the decay of the entering flow
        """

        if self.node_type == "cfd":

            self.cfd_sim.write_vtk(file_path, self.inlet_activity_bq_m3, decay_const)

        return 0

    def calculate_mc_averages(self):
        """
        after reaching the end of the mc sampling, this function calculates
        the mean values
        """

        self.mc_average_activity_bq_m3 = 0
        self.mc_out_activity_bq_m3 = 0
        self.mc_res_time = 0
        self.mc_out_time_cumulated  = 0
        self.tot_activity_bq = 0
        self.linear_velocity_m_s = 0

        #indip_counters = np.count_nonzero(self.mc_counter_vec)
        tot_counts = np.sum(self.mc_counter_vec)

        print ("node", self.node_id)
        #print("indip_counters", indip_counters)
        #print("tot_counts", tot_counts)
        #print ("counter vec", self.mc_counter_vec)
        #print ("average activity vec", self.mc_average_activity_bq_m3_vec)

        if tot_counts == 0:
            return 0

        for i,count in enumerate(self.mc_counter_vec):

            print ("i", i)
            print ("count", count)

            if count != 0:
                print ("average activity", self.mc_average_activity_bq_m3_vec[i]/count)
                self.mc_average_activity_bq_m3 += self.mc_average_activity_bq_m3_vec[i]/count
                self.mc_res_time += self.mc_res_time_vec[i]/count
                self.mc_out_time_cumulated += self.mc_out_time_cumulated_vec[i]/count
                print ("out activity", self.mc_out_activity_bq_m3_vec[i]/count)
                self.mc_out_activity_bq_m3 += self.mc_out_activity_bq_m3_vec[i]/count

        #for i,count in enumerate(self.mc_counter_vec):

        #    if count != 0:
        #        self.mc_average_activity_bq_m3 += self.mc_average_activity_bq_m3_vec[i]*indip_counters
        #        self.mc_res_time += self.mc_res_time_vec[i]*indip_counters
        #        self.mc_out_time_cumulated += self.mc_out_time_cumulated_vec[i]*indip_counters
        #        self.mc_out_activity_bq_m3 += self.mc_out_activity_bq_m3_vec[i]*indip_counters

        #self.mc_average_activity_bq_m3 /= tot_counts
        #self.mc_res_time /= tot_counts
        #self.mc_out_time_cumulated /= tot_counts
        #self.mc_out_activity_bq_m3 /= tot_counts


        #if self.mc_counter != 0:

        #    self.mc_average_activity_bq_m3 += self.mc_average_activity_bq_m3_sum / self.mc_counter
        #    self.mc_res_time = self.mc_res_time_sum / self.mc_counter
        #    self.mc_out_time_cumulated = self.mc_out_time_cumulated_sum / self.mc_counter
        #    self.mc_out_activity_bq_m3 += self.mc_out_activity_bq_m3_sum / self.mc_counter

        vol_m3 = self.volume_cm3 * 1e-6
        self.tot_activity_bq = self.mc_average_activity_bq_m3 * vol_m3

        if self.mc_res_time != 0:
            self.linear_velocity_m_s = self.length_cm / (self.mc_res_time*100)

        return 0


    def increase_mc_counter(self,idx=0):
        """
        this function increases the counter of the node
        """

        self.mc_counter += 1
        self.mc_counter_vec[idx] += 1


        return 0

    def reset_mc_counters(self, n=1):
        """
        this function resets the counters of the node
        """
        self.mc_error = 0


        self.mc_average_activity_bq_m3 = 0
        self.mc_average_activity_bq_m3_sum = 0
        self.mc_average_activity_bq_m3_vec = np.zeros(n)


        self.mc_average_activity_bq_m3_sq_sum = 0
        self.mc_average_activity_bq_m3_sq_vec = np.zeros(n)

        self.mc_counter = 0
        self.mc_counter_vec = np.zeros(n)

        self.mc_res_time = 0
        self.mc_res_time_sum = 0
        self.mc_res_time_vec = np.zeros(n)

        self.mc_out_time_cumulated = 0
        self.mc_out_time_cumulated_sum = 0
        self.mc_out_time_cumulated_vec = np.zeros(n)

        self.mc_out_activity_bq_m3 = 0
        self.mc_out_activity_bq_m3_sum = 0
        self.mc_out_activity_bq_m3_vec = np.zeros(n)


        return 0
