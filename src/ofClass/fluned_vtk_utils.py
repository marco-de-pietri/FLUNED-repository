import pyvista as pv
from pathlib import Path
import numpy as np
import vtk
import meshio
import os
from collections.abc import Sequence


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


def write_cartesian_vtk(
    out_vtk_path: str | Path,
    dims: tuple[int, int, int],
    spacing: tuple[float, float, float],
    origin: tuple[float, float, float],
    dataset=None,
    dataset_name=None,
):
    """
    Write a Cartesian mesh (vtkImageData).
    - Works with **no** dataset (creates an empty grid).
    - Works with a single array & name.
    - Works with multiple arrays & names (must align 1-to-1).

    Parameters
    ----------
    out_vtk_path : str | Path
        Output VTK file.
    dims : (nx+1, ny+1, nz+1)
        Grid point dimensions.
    spacing : (dx, dy, dz)
        Cell spacing.
    origin : (x0, y0, z0)
        Grid origin.
    dataset : np.ndarray | Sequence[np.ndarray] | None
        Cell-data array(s).
        If None  write an empty grid.
    dataset_name : str | Sequence[str] | None
        Name(s) for the dataset(s).
        Required when `dataset` is provided.
        If None and a single dataset is given  defaults to "data".
    """
    grid = pv.ImageData(dimensions=dims, spacing=spacing, origin=origin)

    if dataset is not None:
        # Normalise to sequences
        is_seq = isinstance(dataset, Sequence) and not isinstance(
            dataset, (np.ndarray, bytes, str)
        )
        data_seq = list(dataset) if is_seq else [dataset]

        if dataset_name is None:
            if len(data_seq) == 1:
                names_seq = ["data"]
            else:
                raise ValueError(
                    "`dataset_name` must be provided for multiple datasets"
                )
        else:
            is_seq_name = isinstance(dataset_name, Sequence) and not isinstance(
                dataset_name, (bytes, str)
            )
            names_seq = list(dataset_name) if is_seq_name else [dataset_name]

        if len(data_seq) != len(names_seq):
            raise ValueError("`dataset` and `dataset_name` must have the same length")

        for arr, name in zip(data_seq, names_seq):
            grid.cell_data[name] = np.asarray(arr)

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


def generate_triangularized_h5m_um_mesh(
    original_vtk_path, openmc_h5m_source_file, scale_factor=100
):
    """
    this function scale an original vtk file, empties it, triangularize it and save it as
    an h5m file
    """
    # define a temporary file
    temp_vtk_path = Path(original_vtk_path).parent / "intermediate_tempfile.vtk"

    # initialization and removal of other arrays
    if original_vtk_path.lower().endswith(".vtu"):
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:  # legacy ASCII/Binary *.vtk
        reader = vtk.vtkUnstructuredGridReader()
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
    reader.SetFileName(original_vtk_path)
    reader.Update()
    mesh = reader.GetOutput()
    mesh.GetPointData().Initialize()  # remove all point data
    mesh.GetCellData().Initialize()  # remove all cell data

    # scale mesh
    sx, sy, sz = scale_factor, scale_factor, scale_factor
    tfm = vtk.vtkTransform()
    tfm.Scale(sx, sy, sz)
    tfilter = vtk.vtkTransformFilter()
    tfilter.SetTransform(tfm)
    tfilter.SetInputData(mesh)
    tfilter.Update()
    mesh_scaled = tfilter.GetOutput()
    tri = vtk.vtkDataSetTriangleFilter()
    tri.SetInputData(mesh_scaled)
    tri.SetTetrahedraOnly(True)
    tri.Update()
    mesh_tet = tri.GetOutput()

    # convert to h5m
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileTypeToBinary()
    writer.SetFileName(str(temp_vtk_path))
    if vtk.VTK_MAJOR_VERSION < 6:
        writer.SetInput(mesh_tet)
    else:
        writer.SetInputData(mesh_tet)
    writer.Write()

    meshio_object = meshio.read(temp_vtk_path)
    meshio.write(openmc_h5m_source_file, meshio_object)

    temp_vtk_path.unlink()

    return


def generate_triangularized_scalar_mesh(
    vtk_in_path, vtk_out_path, scaling, cell_data_array
):
    """
    This function scale and triangularize an unstructured mesh
    It removes all the data except for the cell_data_array
    """

    # scale the simulation vtk from meters to cm and triangularize it
    if vtk_in_path.lower().endswith(".vtu"):
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:  # legacy ASCII/Binary *.vtk
        reader = vtk.vtkUnstructuredGridReader()
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()

    reader.SetFileName(vtk_in_path)
    reader.Update()
    mesh = reader.GetOutput()
    mesh.GetPointData().Initialize()  # remove all point data
    cell_data = mesh.GetCellData()
    for i in reversed(range(cell_data.GetNumberOfArrays())):  # iterate safely
        if cell_data.GetArrayName(i) != cell_data_array:
            cell_data.RemoveArray(i)
    sx, sy, sz = scaling, scaling, scaling

    tfm = vtk.vtkTransform()
    tfm.Scale(sx, sy, sz)

    tfilter = vtk.vtkTransformFilter()
    tfilter.SetTransform(tfm)
    tfilter.SetInputData(mesh)
    tfilter.Update()
    mesh_scaled = tfilter.GetOutput()

    tri = vtk.vtkDataSetTriangleFilter()
    tri.SetInputData(mesh_scaled)
    tri.SetTetrahedraOnly(True)
    tri.Update()
    mesh_tet = tri.GetOutput()

    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileTypeToBinary()
    writer.SetFileName(vtk_out_path)
    if vtk.VTK_MAJOR_VERSION < 6:
        writer.SetInput(mesh_tet)
    else:
        writer.SetInputData(mesh_tet)
    writer.Write()

    return
