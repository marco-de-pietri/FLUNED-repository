import pickle
from pathlib import Path

import h5py
import numpy as np
from lxml import etree as et  # type: ignore[attr-defined]

from fluned.ofoam_class.fluned_mesh_utils import generate_triangularized_h5m_um_mesh
from fluned.ofoam_class.ofoam_base.oFoamBase import is_valid_openfoam_case_directory

from .fluned_case_class import flunedCase


def openmc_fluned_coupling(
    filename: str,
    mat_densities,
    mesh_width,
    mesh_dimension,
    mesh_lower_left,
    all_micro_xs,
    flux_in_each_mesh_voxel,
    *,
    compression: str | None = "gzip",
    level: int = 4,
) -> None:
    """
    This function can be imported in a openmc python script to export activation data required for coupled openmc/fluned simulations
    """

    micro_xs_data = np.stack([np.squeeze(a.data, axis=-1) for a in all_micro_xs])
    fluxes = [flux[0] for flux in flux_in_each_mesh_voxel]

    input_dict = {}
    input_dict["isotope_density"] = mat_densities
    input_dict["mesh_width"] = mesh_width
    input_dict["mesh_dimension"] = mesh_dimension
    input_dict["mesh_lower_left"] = mesh_lower_left
    input_dict["fluxes"] = fluxes
    input_dict["index_nuc"] = all_micro_xs[0]._index_nuc
    input_dict["index_rx"] = all_micro_xs[0]._index_rx

    with h5py.File(filename, "w") as f:
        # ── store Python objects ──────────────────────────────────────────────
        for key, value in input_dict.items():
            blob = pickle.dumps(value, protocol=pickle.HIGHEST_PROTOCOL)
            arr = np.frombuffer(blob, dtype=np.uint8)
            f.create_dataset(
                key,
                data=arr,
                dtype=np.uint8,
                compression=compression,
                compression_opts=level,
            )

        # ── store a large array raw ────────────────────────────────
        if micro_xs_data is not None:
            array = micro_xs_data
            ds_name = "data"
            chunk_shape = tuple(min(s, 1024) for s in array.shape)  # simple chunking
            f.create_dataset(
                ds_name,
                data=array,
                dtype=array.dtype,
                chunks=chunk_shape,
                compression=compression,
                compression_opts=level,
            )

    return


class flunedCaseMultiIsotopes:
    """
    this class manages multiple fluned simulations
    """

    def __init__(self, multi_cases_path):
        """
        the init function only gathers the possible folders where fluned simulations are present
        """

        self.container_path = Path(multi_cases_path)
        self.fluned_cases = self.get_fluned_cases()

    def get_fluned_cases(self):
        root = self.container_path
        return [
            flunedCase(fluned_path=p)
            for p in root.iterdir()
            if p.is_dir() and is_valid_openfoam_case_directory(p)
        ]

    def post_process_cases(self):
        """
        post process each fluned case
        """

        for sim in self.fluned_cases:
            sim.parse_fluned_simulation()

        return

    def generate_openmc_um_source(self, mesh_id=100):
        """
        triangularize each unstructured mesh results and generate
        the commands required to import the source model
        """

        xml_source_basename = "source_mesh_file_"
        h5m_basename = "geometry.h5m"
        h5m_file_path = self.container_path / h5m_basename

        commands_file_path = self.container_path / "openmc_source_commands.txt"

        xml_mesh_file_name = "um_mesh_file.xml"
        xml_mesh_file_path = self.container_path / xml_mesh_file_name

        generate_triangularized_h5m_um_mesh(
            self.fluned_cases[0].fluned_simulation.vtk_file_path, h5m_file_path
        )

        for case in self.fluned_cases:
            vtk_file_path_temp = self.container_path / "triangularized.vtk"
            case.fluned_simulation.compute_triangularized_emission_rates(
                vtk_file_path_temp,
                scaling_factor=100.0,
                cell_data_array="T",
                save_tri_mesh_vtk=False,
            )

        # check that all the triangularized scalars have the same length

        if (
            len(
                set(
                    [
                        len(case.fluned_simulation.tri_mesh_emission_rates)
                        for case in self.fluned_cases
                    ]
                )
            )
            > 1
        ):
            print(
                [
                    len(case.fluned_simulation.tri_mesh_emission_rates)
                    for case in self.fluned_cases
                ]
            )
            raise ValueError(
                "ERROR not all the triangularized isotope concentrations are the same"
            )

        #
        # generate xml file with the geometry mesh reference
        #

        # create root element
        root = et.Element("mesh_geo_source")

        # create mesh element with id attribute
        mesh = et.SubElement(
            root,
            "mesh",
            id=str(mesh_id),
            name="source_mesh",
            type="unstructured",
            library="moab",
        )

        # add child elements with text content
        filename = et.SubElement(mesh, "filename")
        filename.text = h5m_basename

        # write to file with xml declaration
        tree = et.ElementTree(root)
        tree.write(
            xml_mesh_file_path,
            encoding="utf-8",
            pretty_print=True,
            xml_declaration=True,
        )

        #
        # generate xml file for each source
        #

        isotope_vec = []

        for case in self.fluned_cases:
            isotope = case.fluned_simulation.isotope.lower()

            isotope_vec.append(isotope)

            case_xml_file = xml_source_basename + isotope + ".xml"
            case_xml_file_path = self.container_path / case_xml_file

            case_root = et.Element("mesh_source")

            # create sublement with the mesh source
            source_mesh = et.SubElement(
                case_root,
                "source",
                type="independent",
                particle=case.fluned_simulation.particle_type,
                strength=f"{case.fluned_simulation.total_isotope_emission_rate:.6e}",
            )

            space = et.SubElement(
                source_mesh,
                "space",
                type="mesh",
                mesh_id=str(mesh_id),
                volume_normalized="False",
            )
            strengths = et.SubElement(space, "strengths")
            strengths.text = " ".join(
                map(
                    str,
                    [val for val in case.fluned_simulation.tri_mesh_emission_rates],
                )
            )  # adjust that the decay rate has been calculated after scaling the vtk file

            ebins_temp = [
                e * 1e6 for e in case.fluned_simulation.e_lines
            ]  # convert from MeV to eV
            energy_parameters = [*ebins_temp, *case.fluned_simulation.p_lines]
            energy = et.SubElement(source_mesh, "energy", type="discrete")
            params = et.SubElement(energy, "parameters")
            params.text = " ".join(map(str, energy_parameters))

            # write to file with xml declaration
            case_tree = et.ElementTree(case_root)
            case_tree.write(
                case_xml_file_path,
                encoding="utf-8",
                pretty_print=True,
                xml_declaration=True,
            )

        with open(commands_file_path, "w") as f:
            f.write("from lxml import etree\n")
            f.write("parser = etree.XMLParser(huge_tree=True)\n")
            f.write(
                'source_root_mesh = etree.parse("{}", parser=parser).getroot()\n'.format(
                    xml_mesh_file_name
                )
            )
            f.write("mesh_element = source_root_mesh.find('mesh')\n")
            f.write(
                "mesh_geo = openmc.UnstructuredMesh.from_xml_element(mesh_element)\n"
            )

            f.write("sources=[]\n")
            source_string = '"' + '","'.join(isotope_vec) + '"'
            f.write(f"isotope_sources = [{source_string}]\n")
            f.write("for isotope in isotope_sources:\n")
            f.write(
                '    xml_import_name = "{}" + isotope + ".xml"\n'.format(
                    xml_source_basename
                )
            )

            f.write(
                "    source_root = etree.parse(xml_import_name, parser=parser).getroot()\n"
            )

            f.write("    source_element = source_root.find('source')\n")
            f.write(
                "    source = openmc.IndependentSource.from_xml_element(source_element, {100:mesh_geo})\n"
            )
            f.write("    sources.append(source)\n")

        return
