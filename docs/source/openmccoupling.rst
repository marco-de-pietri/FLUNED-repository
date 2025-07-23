Multiple Isotopes workflow
--------------------------

FLUNED incorporates the possilibity of using activation data for an arbitrary irradiated flowing material, calculated in openMC, and use it in FLUNED simulations.
This requires to perform an activation calculation with the openmc.deplete module over a cartesian mesh tally covering the flowing material spatial domain. 
The user is responsible for checking the spatial correspondence between the CFD mesh domain and the relevant OpenMC cells. 
However, exact spatial correspondence is not required as reaction rates are computed by combining the neutron fluxes in each mesh voxel with the parent radioisotopes atom densities as provided by the OpenMC material, regardless of the cell positions.
This type of simulation requires the accessibility of the openmc modules in the python environment where the fluned pre- and post-processors are installed.



With this data, the fluned pre-processor generates a FLUNED case for each radioisotope emitting photons that are generated in the flowing fluid.
The calculation assumptions used in single-isotope calculations still hold: only one-step reactions are considered, the parent isotope concentration is kept constant, and the radioisotope concentration is treated as a passive scalar.
The relevant radioisotopes are determined using the chain file used by the openmc.deplete module.

The transfer of data from an openmc simulation to fluned is done by saving the relevant data into an H5 file whose path is passed as the **ACTIVATION_FILE** parameter in the fluned pre-processor input.
The data is packaged by the helper function `openmc_fluned_coupling defined` in the source file: ``FLUNED-Repository/src/fluned_case_multi_isotopes_class.py``.
This function must be imported and called in the openmc input script file after running the function `openmc.deplete.get_microxs_and_flux` function of the deplete module.
An example of this process is shown in the code block below.



code::

    def openmc_fluned_coupling():
        """
        This function can be imported in a openmc python script to run openmc simulations
        """

        return



After running the FLUNED solver for the  various radioisotopes a specific post-processor for the multiple-isotope workflow.


