Multiple Isotopes workflow
--------------------------

FLUNED incorporates the possilibity of using activation data for an arbitrary irradiated flowing fluid, calculated in openMC, and use it in FLUNED simulations.
This requires to perform an activation calculation with the openmc.deplete module over a cartesian mesh tally covering the flowing material spatial domain. 
The user is responsible for checking the spatial correspondence between the CFD mesh domain and the relevant OpenMC cells. 
However, exact spatial correspondence is not required as reaction rates are computed by combining the neutron fluxes in each mesh voxel with the parent radioisotopes atom densities as provided by the OpenMC material, regardless of the cell positions.
This type of simulation requires the accessibility of the openmc modules in the python environment where the fluned pre- and post-processors are installed.

With this data, the fluned pre-processor generates a FLUNED case for each radioisotope emitting photons that are generated in the flowing fluid.
The calculation assumptions used in single-isotope calculations still hold: only one-step reactions are considered, the parent isotope concentration is kept constant, and the radioisotope concentration is treated as a passive scalar.
The relevant radioisotopes are determined using the chain file used by the openmc.deplete module.

The transfer of data from an openmc simulation to fluned is done by saving the relevant data into an H5 file whose path is passed as the **ACTIVATION_FILE** parameter in the fluned pre-processor input.
The data is packaged by the helper function `openmc_fluned_coupling defined` in the source file: ``FLUNED-Repository/src/fluned_case_multi_isotopes_class.py``.
This function can be copied or imported and called in the openmc input script file after running the function `openmc.deplete.get_microxs_and_flux` function of the deplete module.
An example of this process is shown in the code block below.



.. code-block:: python

    import openmc
    import openmc.deplete
    from pathlib import Path
    from openmc_fusion_benchmarks.neutron_sources import fng_source
    from fluned.fluned_case_multi_isotopes_class import openmc_fluned_coupling


    # LOAD GEOMETRY AND MATERIALS

    model_folder = Path("./01_MODEL")

    materials = openmc.Materials.from_xml(model_folder / "materials.xml")


    # assign name to flibe material
    for mat in materials:
        mat.depletable = True
        nuclides = mat.get_nuclides()
        flibe_list = ["Li6", "Li7", "Be9", "F19"] # FLIBE is used for the example
        if all(nuclide in flibe_list for nuclide in nuclides):
            mat.name = "flibe"
            isotope_list = mat.get_nuclides()
            densities = mat.get_nuclide_atom_densities()


    geometry = openmc.Geometry.from_xml(
        path=model_folder / "geometry.xml", materials=materials
    )

    # IMPORT FNG SOURCE

    fng_center = (0, 0, 0.0005)
    fng_uvw = (0.0, 1.0, 0)
    strength = 1.0

    source = fng_source(center=fng_center, reference_uvw=fng_uvw, beam_energy=230)

    # the original strenghts as imported sum to 1
    for s in source:
        s.strength = s.strength * strength



    mesh = openmc.RegularMesh().from_domain(geometry.get_all_cells()[200])


    flux_mesh_size_dm = 0.05
    x_width, y_width, z_width = mesh.width
    mesh.dimension = [
        max(1, int(x_width / flux_mesh_size_dm)),
        max(1, int(y_width / flux_mesh_size_dm)),
        max(1, int(z_width / flux_mesh_size_dm)),
    ]


    mesh_filter = openmc.MeshFilter(mesh)


    # SETTINGS

    settings = openmc.Settings(run_mode="fixed source")
    settings.batches = 10
    settings.inactive = 0
    settings.run_mode = "fixed source"
    settings.source = source
    settings.particles = int(1e5)
    settings.cutoff = {"energy_neutron": 10}
    settings.output = {"tallies": False}

    settings.export_to_xml(path=model_folder / "settings.xml")

    # MODEL GENERATION

    model = openmc.Model()
    model.materials = materials
    model.geometry = geometry
    model.settings = settings


    # RUN

    flux_in_each_mesh_voxel, all_micro_xs = openmc.deplete.get_microxs_and_flux(
        model=model,
        domains=mesh,
        energies=[0, 30e6],
        nuclides=isotope_list,
    )


    openmc_fluned_coupling(
        "fluned_input_file.h5",
        densities,
        mesh.width,
        mesh.dimension,
        mesh.lower_left,
        strength,
        all_micro_xs,
        flux_in_each_mesh_voxel,
    )



After running the FLUNED solver for each radioisotope, a specific post-processor for the multiple-isotope workflow is provided. 
This post-processor can be launched by simply running in the CFD folder containing the various FLUNED cases the following command:

   .. code-block:: bash

      fluned-post-multi-isotopes

This post-processor crawls the folders with the completed FLUNED cases and generates a computational model of the radiation source of the activated fluid.
The geometry of the computational model is the tetrahedralized Unstructured Mesh of the CFD and FLUNED simulations. 
The computational model is composed by a series of xml files, one for each radioisotope and one for the mesh geometry. 
These files can be imported and combined in a OpenMC simulations to define the openmc.settings.source object.
The code to import these is generated automatically by the post-processor and it is saved in the file ``openmc_source_commands.txt``.


