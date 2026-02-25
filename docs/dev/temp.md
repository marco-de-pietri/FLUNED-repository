temporary document to track refactor

foam_case class

- remove default input in the init. use handle input module
- reorganize generate input

     ├╴flunedCase [23, 7]
     ├╴__init__ [30, 9]
     ├╴generate_case [46, 9] refactor

     ├╴generate_fluned_files [242, 9] ok
     ├╴parse_fluned_simulation [292, 9] ok
     ├╴initialize_cfd_class [305, 9] ok

     ├╴create_case_folders_fluent [316, 9] - replaceable with a couple of lines

     ├╴convert_fluent_to_openfoam [330, 9] ok
     ├╴parse_fluent_simulation [339, 9] ok

     ├╴copy_last_phi [352, 9] - definitely simplifiable, maybe can be moved to the ofoamFluned     class?
     ├╴copy_last_u [374, 9]
     ├╴copy_last_nut [386, 9]
     ├╴copy_poly_mesh [400, 9]

     ├╴launch_solver [415, 9] - moveable to the fluned tool launcher

     ├╴launch_check_mode [430, 9] ok

     ├╴generate_cartesian_radiation_source_model [440, 9] okish - could be moved to the ofoamFuned?

     ├╴write_results [456, 9] ok but ugly
     ├╴write_summary_transient [468, 9]
     └╴write_summary_steady [581, 9]



ofoam class


     ├╴formatValues [35, 5]
     ├╴check_float [58, 5]

     ├╴SimulationOF [69, 7]
     ├╴__init__ [80, 9]


     ├╴create_case_folder [104, 9]
     ├╴is_valid_openfoam_case_directory [124, 9]

     ├╴post_process_openfoam_simulation [146, 9] post - both
     ├╴get_number_internal_cells [163, 9] post - both

     ├╴post_process_fluned_simulation [185, 9] - post fluned

     ├╴get_vtk_file [250, 9] - post fluned - vtk

     ├╴get_patches [288, 9] post - both

     ├╴get_reduction_rate [302, 9] - post fluned - could be rationalized but maybe not worth it
     ├╴get_outlet_rr_conc_atoms_m3 [319, 9]
     ├╴get_inlet_td_conc_atoms_m3 [336, 9]
     ├╴get_outlet_t_conc_atoms_m3 [353, 9]
     ├╴get_total_inlet_t_atoms [370, 9]
     ├╴get_total_outlet_t_atoms [384, 9]
     ├╴get_normalized_average_td [398, 9]

     ├╴parse_boundary_phi_files [411, 9] both - optimizable

     ├╴get_last_time [504, 9] - both pre/post replaceable with libfoam
     ├╴get_time_folders [515, 9] both pre/post replaceable with libfoam

     ├╴get_volumetric_flag [538, 9] both pre/post replaceable with libfoam
     ├╴get_density [557, 9] both pre/post replaceable with libfoam
     ├╴get_phi_dimensions [570, 9] both pre/post replaceable with libfoam

     ├╴query_dimensions [582, 9] both pre/post  - block of functions that allow the conversion of  huge files (not replaceable at the moment)
     ├╴query_of_single_value [600, 9] both pre/post 
     ├╴tokenizer [626, 9] both pre/post 
     ├╴token_classifier [648, 9] both pre/post 
     ├╴convert_phi_to_volumetric [779, 9]
     ├╴write_of_tokens [793, 9]
     ├╴phi_tokens_to_volumetric [817, 9]

     ├╴scale_mesh_results [900, 9] post - fluned important

     ├╴calculate_cartesian_sampling_coordinates [934, 9] fluned post - generation of cartesian mesh
     ├╴sample_source_to_cartesian_mesh [1006, 9] fluned post generation of cartesian mesh
     ├╴write_sampled_cartesian_source_vtk [2403, 9] fluned post - this part should be modified to avoid storing so many class attributes relative to this sampling.

     ├╴generate_system_files [1060, 9] fluned - pre prime to refactor with foamlibb
     ├╴generate_constant_file [1524, 9] fluned - pre prime to refactor  with faomlib
     ├╴generate_zero_t [1572, 9]
     ├╴generate_zero_ta [1628, 9]
     ├╴generate_zero_tr [1684, 9]
     ├╴generate_zero_td [1742, 9]

     ├╴read_volumes [1798, 9] both pre/post - should be totally replaced by foamlib as they are not streamed but just read and stored in memory
     ├╴read_centroids [1849, 9]
     ├╴read_velocities [1910, 9]
     ├╴read_grad_t [1938, 9]

     ├╴assign_activation_rates [1964, 9] pre fluned - important

     ├╴generate_source_file [2056, 9] pre fluned - should be totally replaced by foamlib as they are not streamed but just read and stored in memory
     ├╴generate_tr_source_file [2150, 9]

     ├╴assign_isotope_data [2215, 9] - totally need to reevaluate as I am already doing this in the
     isotope class

     ├╴parse_constants_file [2278, 9] - replaceable with foamlib

     ├╴get_time_treatment [2329, 9]
     ├╴read_t [2361, 9]


     ├╴write_cdgs [2436, 9] - fluned post at the moment no modifications
     ├╴write_openmc_sm_source [2508, 9] fluned - post
     ├╴write_openmc_um_source [2607, 9] fluned - post
     └╴compute_triangularized_emission_rates [2709, 9] 
