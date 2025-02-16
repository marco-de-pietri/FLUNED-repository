"""
this module is able to write a vtk file of the pipe circuit
"""
import sys
import math
import numpy as np
import pyvista as pv
import vtk

from .utils import extend_pnt

def orthogonal_vector(axis):
    """
    function to calculate an orthogonal vector to a given vector
    """

    if axis[0]!= 0:
        y_2 = 1
        z_2 = 1
        x_2 = (-axis[1]*y_2-axis[2]*z_2)/axis[0]
    elif axis[1]!= 0:
        x_2 = 1
        z_2 = 1
        y_2 = (-axis[0]*x_2-axis[2]*z_2)/axis[1]
    elif axis[2]!= 0:
        x_2 = 1
        y_2 = 1
        z_2 = (-axis[0]*x_2-axis[1]*y_2)/axis[2]
    else:
        print("wrong axis error")
        sys.exit()

    norm = pow(sum([pow(val,2) for val in [x_2,y_2,z_2]]),0.5)


    unity_norm = [val/norm for val in [x_2,y_2,z_2]]

    #print (axis)
    #print (unityNorm)

    return unity_norm


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a_1 = math.cos(theta / 2.0)
    b_1, c_1, d_1 = -axis * math.sin(theta / 2.0)
    a_a, b_b, c_c, d_d = a_1 * a_1, b_1 * b_1, c_1 * c_1, d_1 * d_1
    b_c, a_d, a_c, a_b, b_d, c_d = ( b_1 * c_1,
                                     a_1 * d_1,
                                     a_1 * c_1,
                                     a_1 * b_1,
                                     b_1 * d_1,
                                     c_1 * d_1)
    res = np.array([[a_a + b_b - c_c - d_d, 2 * (b_c + a_d), 2 * (b_d - a_c)],
                     [2 * (b_c - a_d), a_a + c_c - b_b - d_d, 2 * (c_d + a_b)],
                     [2 * (b_d + a_c), 2 * (c_d - a_b), a_a + d_d - b_b - c_c]])
    return res


def write_vtk_files(nodes,full_path, decay_const, steady_state_flag):
    """
    this function writes the vtk file of the pipe circuit for visualization
    """

    #pipes = {key:value for key,value in nodes.items() if value.node_type in ["pipe","tank-cyl"]}

    polyhedron_nodes = {key:value for key,value in nodes.items() if value.node_type in ["pipe","tank-cyl","stl"]}

    cfd_nodes = {key:value for key,value in nodes.items() if value.node_type in ["cfd"]}

    #write_vtk_pipes(pipes,full_path, steady_state_flag)

    write_vtk_pipes_polyhedron(polyhedron_nodes, full_path, steady_state_flag )

    write_vtk_cfd_nodes(cfd_nodes,full_path, steady_state_flag, decay_const)

    return 0

def write_vtk_cfd_nodes(cfd_nodes,full_path, steady_state_flag, decay_const):
    """
    this function writes the vtk file of the cfd nodes for visualization
    """

    if not steady_state_flag and len(cfd_nodes) > 0:
        print ("vtk of the cfd is not produced in transient mode")
        return



    for node in cfd_nodes.values():

        node_path = full_path + "_cfd_node_" + str(node.node_id) + ".vtk"

        node.write_cfd_vtk(node_path, decay_const)






    return


def write_vtk_pipes(pipes,full_path, steady_state_flag):
    """
    this function writes the vtk file of the pipe nodes for visualization
    """

    full_path = full_path + "_pipes.vtk"

    n_pipes = len(pipes)
    try:
        out = open(full_path,'w',encoding="utf8", errors='ignore')
    except IOError:
        print("couldn't write vtk file")
        sys.exit()
    with out:

        header = """# vtk DataFile Version 2.0
DGS printing
ASCII
DATASET UNSTRUCTURED_GRID\r\n"""
        out.write(header)
        out.write(f"POINTS  {(n_pipes*24):d}  floats\r\n")

        for pipe in pipes.values():
            pipe_axis = [pipe.axis_x, pipe.axis_y, pipe.axis_z]
            pipe_axis_1 = [pipe.axis_x_1, pipe.axis_y_1, pipe.axis_z_1]
            radius_m = pipe.radius_cm/100
            length_m = pipe.length_cm/100
            pipe_origin_m = [pipe.origin_x_cm/100,
                             pipe.origin_y_cm/100,
                             pipe.origin_z_cm/100]
            pipe_end_m = extend_pnt(pipe_origin_m,
                                           pipe_axis,
                                           length_m)
            pipe_mid_m = extend_pnt(pipe_origin_m,
                                           pipe_axis,
                                           length_m/2)

            vector0   = pipe_axis_1

            pi_val = 3.141592653589793

            vector45  = np.dot(rotation_matrix(pipe_axis,pi_val/4), vector0)
            vector90  = np.dot(rotation_matrix(pipe_axis,pi_val/2), vector0)
            vector135 = np.dot(rotation_matrix(pipe_axis,pi_val*3/4),vector0)
            vector180 = np.dot(rotation_matrix(pipe_axis,pi_val), vector0)
            vector225 = np.dot(rotation_matrix(pipe_axis,pi_val*5/4), vector0)
            vector270 = np.dot(rotation_matrix(pipe_axis,pi_val + pi_val/2),
                               vector0)
            vector315 = np.dot(rotation_matrix(pipe_axis,pi_val*7/4), vector0)

            # reference :   https://vtk.org/doc/release/5.2/html/a00103.html

            point00 = extend_pnt(pipe_origin_m,vector0,  radius_m)
            point01 = extend_pnt(pipe_origin_m,vector90, radius_m)
            point02 = extend_pnt(pipe_origin_m,vector180,radius_m)
            point03 = extend_pnt(pipe_origin_m,vector270,radius_m)
            point04 = extend_pnt(pipe_end_m,vector0,  radius_m)
            point05 = extend_pnt(pipe_end_m,vector90, radius_m)
            point06 = extend_pnt(pipe_end_m,vector180,radius_m)
            point07 = extend_pnt(pipe_end_m,vector270,radius_m)
            point08 = extend_pnt(pipe_origin_m,vector45,radius_m)
            point09 = extend_pnt(pipe_origin_m,vector135,radius_m)
            point10 = extend_pnt(pipe_origin_m,vector225,radius_m)
            point11 = extend_pnt(pipe_origin_m,vector315,radius_m)
            point12 = extend_pnt(pipe_end_m,vector45, radius_m)
            point13 = extend_pnt(pipe_end_m,vector135,radius_m)
            point14 = extend_pnt(pipe_end_m,vector225,radius_m)
            point15 = extend_pnt(pipe_end_m,vector315,radius_m)
            point16 = extend_pnt(pipe_mid_m,vector0,  radius_m)
            point17 = extend_pnt(pipe_mid_m,vector90, radius_m)
            point18 = extend_pnt(pipe_mid_m,vector180,radius_m)
            point19 = extend_pnt(pipe_mid_m,vector270,radius_m)
            point20 = extend_pnt(pipe_mid_m,vector315,radius_m)
            point21 = extend_pnt(pipe_mid_m,vector135,radius_m)
            point22 = extend_pnt(pipe_mid_m,vector45, radius_m)
            point23 = extend_pnt(pipe_mid_m,vector225,radius_m)
            pipe_points    = [point00,point01,point02,point03,point04,point05,
                              point06,point07,point08,point09,point10,point11,
                              point12,point13,point14,point15,point16,point17,
                              point18,point19,point20,point21,point22,point23]
            i = 0
            write_string = ""
            for coord in pipe_points:
                if (i%3 == 0 ) and (i > 0):
                    out.write(write_string + '\r\n')
                    write_string = (f"{coord[0]:>.5f} "
                                    f"{coord[1]:>.5f} "
                                    f"{coord[2]:>.5f} ")
                    i = 1
                else:
                    new_string = (f"{coord[0]:>.5f} "
                                  f"{coord[1]:>.5f} "
                                  f"{coord[2]:>.5f} ")
                    write_string += new_string
                    i += 1
            out.write(write_string + '\r\n')

        out.write(f"CELLS  {n_pipes:d}   {(n_pipes*25):d}\r\n")
        for i,pipe in enumerate(pipes):
            pnt = 24*i
            write_string = (f"24 {(pnt+ 0):d} {(pnt+ 1):d} {(pnt+ 2):d} "
                               f"{(pnt+ 3):d} {(pnt+ 4):d} {(pnt+ 5):d} "
                               f"{(pnt+ 6):d} {(pnt+ 7):d} {(pnt+ 8):d} "
                               f"{(pnt+ 9):d} {(pnt+10):d} {(pnt+11):d} "
                               f"{(pnt+12):d} {(pnt+13):d} {(pnt+14):d} "
                               f"{(pnt+15):d} {(pnt+16):d} {(pnt+17):d} "
                               f"{(pnt+18):d} {(pnt+19):d} {(pnt+20):d} "
                               f"{(pnt+21):d} {(pnt+22):d} {(pnt+23):d}\r\n")
            out.write(write_string)

        out.write(f"CELL_TYPES  {n_pipes:d} \r\n")
        for i,pipe in enumerate(pipes):
            out.write("33\r\n")

        out.write(f"CELL_DATA  {n_pipes:d} \r\n")

        out.write("SCALARS node_id int 1 \r\n")
        out.write("LOOKUP_TABLE default\r\n")
        for pipe in pipes.values():
            val = pipe.node_id
            out.write(f"{val:d}\r\n")

        out.write("SCALARS inlet_activity_bq_m3 double 1 \r\n")
        out.write("LOOKUP_TABLE default\r\n")
        for pipe in pipes.values():
            activity = pipe.inlet_activity_bq_m3
            out.write(f"{activity:.5e}\r\n")

        out.write("SCALARS average_vol_activity_bq_m3 double 1 \r\n")
        out.write("LOOKUP_TABLE default\r\n")
        for pipe in pipes.values():
            activity = pipe.average_activity_bq_m3
            out.write(f"{activity:.5e}\r\n")

        out.write("SCALARS total_activity_bq double 1 \r\n")
        out.write("LOOKUP_TABLE default\r\n")
        for pipe in pipes.values():
            tot_activity = pipe.tot_activity_bq
            out.write(f"{tot_activity:.5e}\r\n")

        out.write("SCALARS mc_error double 1 \r\n")
        out.write("LOOKUP_TABLE default\r\n")
        for pipe in pipes.values():
            error = pipe.mc_error
            out.write(f"{error:.5f}\r\n")

        if steady_state_flag is True:

            out.write("SCALARS reynolds double 1 \r\n")
            out.write("LOOKUP_TABLE default\r\n")
            for pipe in pipes.values():
                reynolds = pipe.reynolds
                out.write(f"{reynolds:.5f}\r\n")

            out.write("SCALARS res_time_s double 1 \r\n")
            out.write("LOOKUP_TABLE default\r\n")
            for pipe in pipes.values():
                time = pipe.mc_res_time
                out.write(f"{time:.5f}\r\n")

            out.write("SCALARS mass_flow_kg_s double 1 \r\n")
            out.write("LOOKUP_TABLE default\r\n")
            for pipe in pipes.values():
                mflow = pipe.mass_flow
                out.write(f"{mflow:.5e}\r\n")

            out.write("SCALARS linear_velocity_m_s double 1 \r\n")
            out.write("LOOKUP_TABLE default\r\n")
            for pipe in pipes.values():
                mflow = pipe.linear_velocity_m_s
                out.write(f"{mflow:.5e}\r\n")

    return 0




def create_cylinder_polyhedron(pipe, resolution=24):
    """
    Create a cylinder as a single vtkPolyhedron cell.
    The cylinder is built from two rings of points (bottom and top)
    and the cell is defined by resolution lateral quadrilaterals plus
    one bottom and one top polygon face.
    """
    # Normalize pipe axis and convert dimensions (cm -> m)
    axis = np.array([pipe.axis_x, pipe.axis_y, pipe.axis_z], dtype=float)
    norm = np.linalg.norm(axis)
    if norm == 0:
        raise ValueError("Pipe axis is zero-length.")
    axis /= norm
    radius = pipe.radius_cm / 100.0
    height = pipe.length_cm / 100.0
    origin = np.array([pipe.origin_x_cm, pipe.origin_y_cm, pipe.origin_z_cm], dtype=float) / 100.0

    # Build an orthonormal basis for the circle.
    arbitrary = np.array([1, 0, 0], dtype=float)
    if np.allclose(axis, arbitrary):
        arbitrary = np.array([0, 1, 0], dtype=float)
    v = np.cross(axis, arbitrary)
    v /= np.linalg.norm(v)
    w = np.cross(axis, v)

    # Generate the bottom ring points
    theta = np.linspace(0, 2*np.pi, resolution, endpoint=False)
    bottom_ring = np.array([origin + radius * (np.cos(t) * v + np.sin(t) * w)
                            for t in theta])
    # The top ring is translated along the axis
    top_ring = bottom_ring + axis * height

    # Combine points: bottom ring indices 0...resolution-1, top ring indices resolution...2*resolution-1
    pts = np.vstack((bottom_ring, top_ring))

    # Define faces for the polyhedron:
    faces = []
    # Lateral faces (each quad uses one adjacent pair of bottom and top ring points)
    for i in range(resolution):
        i_next = (i + 1) % resolution
        faces.append([i, i_next, i_next + resolution, i + resolution])
    # Bottom cap (using bottom ring in order)
    faces.append(list(range(0, resolution)))
    # Top cap (using top ring in reverse order to maintain outward normals)
    faces.append(list(range(resolution, 2 * resolution))[::-1])

    nfaces = len(faces)  # should equal resolution + 2

    # Create vtkPoints from pts
    vtk_pts = vtk.vtkPoints()
    for p in pts:
        vtk_pts.InsertNextPoint(p.tolist())

    # Build the connectivity list for the vtkPolyhedron.
    # The required format is: [nfaces, n0, id0, id1, ..., n1, id0, id1, ...]
    id_list = vtk.vtkIdList()
    id_list.InsertNextId(nfaces)
    for face in faces:
        id_list.InsertNextId(len(face))
        for pid in face:
            id_list.InsertNextId(pid)

    # Create an empty unstructured grid and add the polyhedron cell.
    ug = vtk.vtkUnstructuredGrid()
    ug.SetPoints(vtk_pts)
    # vtk.VTK_POLYHEDRON is the cell type for a polyhedron (its integer value is 42)
    ug.InsertNextCell(vtk.VTK_POLYHEDRON, id_list)

    return pv.wrap(ug)

def load_stl_polyhedron(stl_path):
    """
    Load an STL file and return its mesh wrapped as a PyVista object.
    Assumes the STL file represents a single closed polyhedron.
    If the loaded mesh contains more than one cell, all cells are assigned
    the same data later.
    """
    mesh = pv.read(stl_path)
    # Optionally, you might want to clean or merge the mesh here if needed:
    # mesh = mesh.clean()
    return mesh





def write_vtk_pipes_polyhedron(pipes, full_path, steady_state_flag, resolution=24):
    """
    For each pipe node, create a polyhedron cell.
    If pipe.node_type is 'stl', load the geometry from pipe.stl_path.
    Otherwise, build a cylinder polyhedron.
    Then assign the pipe properties as cell data and merge all into one grid,
    which is saved to file.
    """
    full_path = full_path + ".vtk"
    merged = None

    for pipe in pipes.values():
        # Decide which geometry to use based on node_type.
        if hasattr(pipe, "node_type") and pipe.node_type == 'stl':
            try:
                cell_mesh = load_stl_polyhedron(pipe.stl_path)
            except Exception as e:
                print(f"Error loading STL for node {pipe.node_id}: {e}")
                continue
        else:
            cell_mesh = create_cylinder_polyhedron(pipe, resolution=resolution)

        # For consistency, assign the same property value to every cell in this mesh.
        n_cells = cell_mesh.n_cells
        cell_mesh.cell_data["node_id"] = np.full(n_cells, pipe.node_id, dtype=int)
        cell_mesh.cell_data["inlet_vol_activity_bq_m3"] = np.full(n_cells, pipe.inlet_activity_bq_m3, dtype=float)
        cell_mesh.cell_data["average_vol_activity_bq_m3"] = np.full(n_cells, pipe.average_activity_bq_m3, dtype=float)
        cell_mesh.cell_data["total_activity_bq"] = np.full(n_cells, pipe.tot_activity_bq, dtype=float)
        cell_mesh.cell_data["reaction_rate_m3"] = np.full(n_cells, pipe.reaction_rate_m3, dtype=float)
        if steady_state_flag:
            cell_mesh.cell_data["reynolds"] = np.full(n_cells, pipe.reynolds, dtype=float)
            cell_mesh.cell_data["res_time_s"] = np.full(n_cells, pipe.mc_res_time, dtype=float)
            cell_mesh.cell_data["mass_flow_kg_s"] = np.full(n_cells, pipe.mass_flow, dtype=float)
            cell_mesh.cell_data["linear_velocity_m_s"] = np.full(n_cells, pipe.linear_velocity_m_s, dtype=float)

        # Merge this cell into the overall grid.
        if merged is None:
            merged = cell_mesh
        else:
            merged = merged.merge(cell_mesh)

    if merged is None:
        print("No pipes to write.")
        return 0

    # Save the merged unstructured grid containing one polyhedron cell per node.
    merged.save(full_path)
    return 0