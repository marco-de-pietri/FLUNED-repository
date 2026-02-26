from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Literal, Tuple

import meshio
import numpy as np
import pyvista as pv
import vtk

EulerOrder = Literal["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]
CenterMode = Literal["origin", "mesh_center", "bounds_center"]


def generate_external_stl(vtk_path, output_folder):
    """
    this function extracts the external surface of the final vtk file and
    write it in a stl file
    """

    stl_path = Path(output_folder) / "external_surface.stl"

    # read the VTK mesh (handles .vtk, .vtu, .vti, etc.)
    mesh = pv.read(str(vtk_path))

    # extract the exterior faces (drops internal cells)
    surface = mesh.extract_surface(algorithm=None).triangulate()

    # write STL (binary by default)
    surface.save(str(stl_path))

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


def _euler_rotation_matrix(
    angles_deg: tuple[float, float, float],
    order: EulerOrder = "xyz",
) -> np.ndarray:
    """Return 3x3 rotation matrix for extrinsic Euler rotations applied in given order.

    Angles are in degrees.

    Convention here (column vectors): applying rotations sequentially means
    p' = R_order @ p, where for order='xyz': p' = Rz @ Ry @ Rx @ p.
    """
    ax, ay, az = np.deg2rad(angles_deg)

    cx, sx = np.cos(ax), np.sin(ax)
    cy, sy = np.cos(ay), np.sin(ay)
    cz, sz = np.cos(az), np.sin(az)

    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]], dtype=float)
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]], dtype=float)
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]], dtype=float)

    mats = {"x": Rx, "y": Ry, "z": Rz}
    # Apply in order: p' = R_last @ ... @ R_first @ p
    R = np.eye(3, dtype=float)
    for axis in order:
        R = mats[axis] @ R
    return R


def _pick_center(mesh: pv.DataSet, center: CenterMode) -> np.ndarray:
    if center == "origin":
        return np.zeros(3, dtype=float)
    if center == "mesh_center":
        return np.asarray(mesh.center, dtype=float)
    if center == "bounds_center":
        b = mesh.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)
        return np.array(
            [(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, (b[4] + b[5]) / 2], dtype=float
        )


def apply_roto_translation_to_vtk_grid(
    input_path: str | Path,
    output_path: str | Path,
    rotation_deg: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    translation: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    *,
    order: EulerOrder = "xyz",
    center: CenterMode = "origin",
    transform_all_input_vectors: bool = False,
) -> pv.DataSet:
    """Read a VTK file, apply rotation+translation, and write it back.

    Supports:
      - pv.StructuredGrid
      - pv.RectilinearGrid
      - pv.ImageData (UniformGrid)

    Transform:
        p' = R (p - c) + c + t

    Notes:
      - Rotating ImageData/RectilinearGrid requires converting to StructuredGrid first
        (their implicit geometry cannot represent arbitrary rotations). :contentReference[oaicite:2]{index=2}
      - For translation-only, we keep original types: ImageData uses .origin, and
        RectilinearGrid uses .x/.y/.z because points aren't settable. :contentReference[oaicite:3]{index=3}
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    mesh = pv.read(str(input_path))
    t = np.asarray(translation, dtype=float)
    rot = np.asarray(rotation_deg, dtype=float)
    do_rotate = not np.allclose(rot, 0.0)

    # --- Translation-only fast path (keeps original grid class) ---
    if not do_rotate:
        if isinstance(mesh, pv.ImageData):
            out = mesh.copy(deep=True)
            out.origin = tuple(np.asarray(out.origin, dtype=float) + t)
            out.save(str(output_path))
            return out

        if isinstance(mesh, pv.RectilinearGrid):
            out = mesh.copy(deep=True)
            out.x = np.asarray(out.x, dtype=float) + t[0]
            out.y = np.asarray(out.y, dtype=float) + t[1]
            out.z = np.asarray(out.z, dtype=float) + t[2]
            out.save(str(output_path))
            return out

        # StructuredGrid or anything with mutable points
        if hasattr(mesh, "points"):
            out = mesh.copy(deep=True)
            out.points = np.asarray(out.points, dtype=float) + t
            out.save(str(output_path))
            return out

        raise TypeError(
            f"Unsupported dataset type for translation-only: {type(mesh).__name__}"
        )

    # --- Rotation (and translation) path ---
    # Convert rectilinear/uniform grids to StructuredGrid first.
    if isinstance(mesh, (pv.ImageData, pv.RectilinearGrid)):
        mesh = (
            mesh.cast_to_structured_grid()
        )  # supported by both :contentReference[oaicite:4]{index=4}
    elif not isinstance(mesh, pv.StructuredGrid):
        # You can broaden this if you want, but you asked specifically for these 3 grid types.
        raise TypeError(
            f"Rotation requested but dataset is {type(mesh).__name__}; "
            "expected StructuredGrid, RectilinearGrid, or ImageData."
        )

    c = _pick_center(mesh, center)
    R = _euler_rotation_matrix(rotation_deg, order=order)

    # Homogeneous 4x4: p' = R p + (t + c - R c)
    M = np.eye(4, dtype=float)
    M[:3, :3] = R
    M[:3, 3] = t + c - (R @ c)

    out = mesh.transform(
        M,
        transform_all_input_vectors=transform_all_input_vectors,
        inplace=False,
    )
    out.save(str(output_path))
    return


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
        If None write an empty grid.
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
    if original_vtk_path.suffix.lower() == ".vtu":
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
    # the scaling factor is to scale the size (from m to cm), not the values
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
        writer.SetInput(mesh_tet)  # type: ignore[attr-defined]
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
    if vtk_in_path.suffix.lower() == ".vtu":
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
        writer.SetInput(mesh_tet)  # type: ignore[attr-defined]
    else:
        writer.SetInputData(mesh_tet)
    writer.Write()

    return
