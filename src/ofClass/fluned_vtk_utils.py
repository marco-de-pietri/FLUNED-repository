import pyvista as pv
import os
import numpy as np


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


def get_vtk_volumes(file_path: str) -> np.ndarray:
    """
    Read a mesh from the given file and return an array of volumes for each cell.

    Parameters
    ----------
    file_path : str
        Path to the mesh file (VTK, VTP, STL, etc.)

    Returns
    -------
    numpy.ndarray
        1D array of cell volumes.
    """
    # Load the mesh
    mesh = pv.read(file_path)
    # Compute cell sizes (volumes)
    mesh = mesh.compute_cell_sizes(volume=True)
    # Extract the 'Volume' array from cell data
    volumes = mesh.cell_data["Volume"]
    return volumes


def write_cartesian_vtk(out_vtk_path, dims, spacing, origin, dataset, dataset_name):
    """
    This function writes a cartesian mesh
    """
    grid = pv.ImageData(dimensions=dims, spacing=spacing, origin=origin)
    grid.cell_data[dataset_name] = dataset
    grid.save(out_vtk_path)

    return


def get_vtk_celldata_array(file_path: str, array_name: str) -> np.ndarray:
    """
    Read a mesh from the given file and return the specified cell-data array.

    Parameters
    ----------
    file_path : str
        Path to the mesh file (VTK, VTP, STL, etc.).
    array_name : str
        Name of the cell-data array to extract.

    Returns
    -------
    numpy.ndarray
        1D array of the requested cell-data values.

    Raises
    ------
    ValueError
        If the specified array_name is not found in the mesh's cell data.
    """
    mesh = pv.read(file_path)
    if array_name not in mesh.cell_data:
        raise ValueError(
            f"Cell-data array '{array_name}' not found. Available arrays: {list(mesh.cell_data.keys())}"
        )
    return mesh.cell_data[array_name]
