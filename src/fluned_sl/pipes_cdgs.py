"""
code to write the CDGS for cylindrical pipes - author J. Alguacil UNED, M. De Pietri UNED
"""

def write_vec_f(vec,npos):
    """
    function description
    """
    j = -1
    line = ''
    for el_val in vec:
        j = j + 1
        if j == npos:
            line = line+'\n'
            j = 0
        line = line + f"{el_val:12.4f} "
    return line

def write_vec_e(vec,npos):
    """
    function description
    """
    j = -1
    line = ''
    for el_val in vec:
        j = j + 1
        if j == npos:
            line = line+'\n'
            j = 0
        line = line + f"{el_val:12.4e} "
    return line

def get_header(nodes):
    """
    calculate the cdgs header factor
    """
    # Get Source Factor
    nmesh   = 0
    act_tot = 0.0
    for node in nodes.values():
        act = node.tot_activity_bq
        act_tot = act_tot + act
        nmesh = nmesh + 1
    return nmesh,act_tot

def write_source_data(pipe,yld,p_bins):
    """
    function description
    """
    line = "source_data\n"


    # index of the mesh, i_tot (ph/s), vol(cm**3), ncel
    vol_cyl = pipe.volume_cm3

    pipe_emission_rate       = pipe.tot_activity_bq*yld
    stat_error = pipe.mc_error

    line = line + f"1  {pipe_emission_rate:12.8e}  {vol_cyl:12.8e} 1\n"

    # [cell  relative_vol  emission_rate]
    line = line + f"{pipe.mcnp_cell} 1.00000000e+00  {pipe_emission_rate:12.8e} \n"
    gammas = []
    error  = []
    for prob in p_bins:
        gammas.append(prob*pipe_emission_rate)
        if prob == 0.0:
            error.append(0.0)
        else:
            error.append(stat_error)

    line = line + write_vec_e(gammas,6) + "\n"
    line = line + write_vec_e(error,6)  + "\n"
    line = line + "end_source_data\n"

    return line

def write_pipe_mesh_structure(pipe):
    """
    write the cylinder mesh of the pipe
    """
    radius = pipe.radius_cm

    length = pipe.length_cm

    origin = [
                pipe.origin_x_cm,
                pipe.origin_y_cm,
                pipe.origin_z_cm,
             ]

    axs =    [
                pipe.axis_x,
                pipe.axis_y,
                pipe.axis_z,
             ]
    vec_2 =  [
                pipe.axis_x_1,
                pipe.axis_y_1,
                pipe.axis_z_1,
             ]


    # Definition of the boundaries
    theta_vals = [0.0,1.0]
    z_vals     = [0.0,length]
    r_vals     = [0.0,radius]


    # Assumed cyl [r, theta,z]
    line = "mesh_type cyl\n"
    line = line + "mesh_boundaries  2 2 2\n"

    # Origin and transformations
    line = line + f"{origin[0]:12.6e}  {origin[1]:12.6e}  {origin[2]:12.6e}\n"
    line = line + f"{axs[0]:12.6e}  {axs[1]:12.6e}  {axs[2]:12.6e}\n"
    line = line + f"{vec_2[0]:12.6e}  {vec_2[1]:12.6e}  {vec_2[2]:12.6e}\n"

    # Mesh boundaries r,theta,z

    line = line + write_vec_f(r_vals,6) + "\n"
    line = line + write_vec_f(theta_vals,6) + "\n"
    line = line + write_vec_f(z_vals,6)     + "\n"

    return line

def write_general_photon_info(p_yield,e_bins):
    """
    this function writes the energy spectrum of the pipes
    """
    line = f"total_source {p_yield:12.9e} \n"
    line = line + "energy_type bins\n"    # lines are not yet supported
    line = line + f"energy_boundaries {len(e_bins):d}\n"
    line = line + write_vec_ff(e_bins,6) +"\n"
    return line


def format_values_cdgs(vector):
    """
    function to format string for the cdgs file
    """
    max_len = 70
    return_st = ''
    new_line = ''

    for item in vector:
        new_number = f'{item:.7e}'
        if return_st == '' and new_line == '':
            new_line = new_number
            continue

        if len(new_line + ' ' + new_number) > max_len :
            return_st +=new_line + '\n'
            new_line = new_number

        else:
            new_line = new_line + ' ' + new_number

    return_st += new_line + '\n'

    return return_st

def make_pipes_cdgs(pipes,dec_lambda, dec_yield,e_bins,p_bins, full_path):
    """
    function to write the CDGS for cylindrical pipes
    """

    # precision in cm of the cfd sampling - will be implemented as main argument later
    precision_cm = 0.5 # cm


    # Open File
    with open(full_path,"w",encoding="utf8") as out:
    # Get Isotope Features

        # Data of the header
        nmesh,act_tot = get_header(pipes)
        source_factor = act_tot * dec_yield
        line = f"num_meshes {nmesh:d} \nglobal_source {source_factor:9.7e} \n"
        out.write(line)


        # Write Each pipe CDGS
        nmesh = 0
        for pipe in pipes.values():
            if pipe.node_type == 'tank':
                print ("WARNING: tank node not included in the CDGS")
                continue

            if pipe.node_type in  ["pipe", "tank-cyl"]:
                act = pipe.tot_activity_bq
                y_pipe = act*dec_yield

                #Introduction of the common decay gamma source
                nmesh = nmesh + 1
                line  = f"mesh_id {nmesh:d}\n"
                line  = line +f"Cylinder Cell {pipe.node_id:d}\n"
                line  = line +"Cooling_time  0.0\n"
                out.write(line)

                #General Photon Info
                line = write_general_photon_info(y_pipe,e_bins)
                out.write(line)

                # Mesh Structure
                line = write_pipe_mesh_structure(pipe)
                out.write(line)

                # Source Data
                line = write_source_data(pipe,dec_yield,p_bins)
                out.write(line)

            if pipe.node_type == 'cfd':
                act = pipe.tot_activity_bq
                y_pipe = act*dec_yield

                #Introduction of the common decay gamma source
                nmesh = nmesh + 1
                line  = f"mesh_id {nmesh:d}\n"
                line  = line +f"Cfd node {pipe.node_id:d}\n"
                line  = line +"Cooling_time  0.0\n"
                out.write(line)

                #General Photon Info
                line = write_general_photon_info(y_pipe,e_bins)
                out.write(line)

                # Mesh Structure
                #xyz_ints, xyz_nodes, sample_vals = calculate_sampling_coordinates(
                #                    pipe.cfd_sim.vtk_dimensions)
                xyz_nodes, sample_vals = pipe.cfd_sim.sample_scaled_vtk(precision_cm)

                tot_activity_sampled = sum([val['average_vol_activity_bq_m3']*
                                            val['volume_cm3']*1e-6
                                            for val in sample_vals])

                # scaling factor to ensure that the total activity is the same

                scaling_factor = act/tot_activity_sampled

                for val in sample_vals:
                    val['voxel_emission'] = (val['average_vol_activity_bq_m3']*
                                            val['volume_cm3']*1e-6*
                                            scaling_factor*
                                            dec_yield)

                print ("tot activity sampled",tot_activity_sampled)
                print ("tot activity fluned_sl",act)
                print ("scaling factor",scaling_factor)

                nx_nodes = len(xyz_nodes[0])
                ny_nodes = len(xyz_nodes[1])
                nz_nodes = len(xyz_nodes[2])


                out.write("mesh_type rec\n")
                out.write(f"mesh_boundaries {nx_nodes:d} {ny_nodes:d} {nz_nodes:d}\n")

                out.write("0.000000e+00  0.000000e+00  0.000000e+00\n")
                out.write("1.000000e+00  0.000000e+00  0.000000e+00\n")
                out.write("0.000000e+00  1.000000e+00  0.000000e+00\n")
                x_string = format_values_cdgs(xyz_nodes[0])
                out.write(x_string)
                y_string = format_values_cdgs(xyz_nodes[1])
                out.write(y_string)
                z_string = format_values_cdgs(xyz_nodes[2])
                out.write(z_string)
                out.write("source_data\n")


                spec_err_string = format_values_cdgs([0]*(len(e_bins)-1))
                voxel_string_1 = "{:d} {:.5e} {:.5e} 1\n"
                voxel_string_2 = "0 1.0 {:.5e}\n"

                for voxel in sample_vals:
                    if voxel['voxel_emission'] == 0.0:
                        continue
                    # [index of the mesh, emission_tot (ph/s), vol(cm**3), ncel]
                    out.write(voxel_string_1.format(voxel['id'],
                                                 voxel['voxel_emission'],
                                                 voxel['volume_cm3']))

                    # [mcnp cell  relative_vol  emission (ph/s)]
                    out.write(voxel_string_2.format(voxel['voxel_emission']))

                    emitting_spectrum =  ([prob*voxel['voxel_emission']
                                         for prob in  p_bins])
                    spectrum_string=format_values_cdgs(emitting_spectrum)
                    out.write(spectrum_string)
                    out.write(spec_err_string)

                out.write("end_source_data\n")



    return 0


def write_vec_ff(vec,npos):
    """
    line formatting for the energy bins
    """
    j = -1
    line = ''
    for el_val in vec:
        j = j + 1
        if j == npos:
            line = line+'\n'
            j = 0
        line = line + f"{el_val:12.7e} "
    return line
