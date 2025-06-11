import pyvista as pv
import os


def generate_external_stl(vtk_path, output_folder):
    """
    this function extracts the external surface of the final vtk file and
    write it in a stl file
    """

    stl_path = os.path.join(output_folder, "external_surface.stl")

    # read the VTK mesh (handles .vtk, .vtu, .vti, etc.)
    mesh = pv.read(vtk_path)

    # extract the exterior faces (drops internal cells)
    surface = mesh.extract_surface().triangulate()

    # write STL (binary by default)
    surface.save(stl_path)

    return


def sample_vtk(vtk_file, dataset_name, coordinates):
    """
    this function reads the vtk and sample some datasets - official function
    """

    vtk_mesh = pv.read(vtk_file)

    sample_points = pv.PolyData(coordinates)

    sampled_values = sample_points.sample(
        vtk_mesh,
        pass_cell_data=True,
        pass_point_data=True,
        progress_bar=True,
    )
    sampled_dataset = sampled_values.point_data[dataset_name]

    return sampled_dataset


def get_vtk_dimensions(vtk_path):
    """
    this function reads the unstructured mesh and get the cartesian
    boundaries
    """

    mesh = pv.read(vtk_path)

    bounds = mesh.bounds
    volume_m3 = mesh.volume

    return bounds, volume_m3


def write_cartesian_vtk(out_vtk_path, dims, spacing, origin, dataset, dataset_name):
    """
    This function writes a cartesian mesh
    """
    grid = pv.ImageData(dimensions=dims, spacing=spacing, origin=origin)
    grid.cell_data[dataset_name] = dataset
    grid.save(out_vtk_path)

    return
