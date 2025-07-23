================
Input File Guide
================

This page documents all keywords that can appear in a `FLUNED` input, their
purpose, default behaviour and the range of values they accept.

To start, a generic template can be created by running `fluned -t`.

Every directive is case‑insensitive, except for the files and folder paths, and must appear on its own line.  
Comment text after a `#` is ignored.

The fluned pre-processor input file can be structured to initiate the generation of more than one case at the time by including multiple sets of input keywords, one for each case.
The keywords set must be separated by at least one empty line and the CASE parameters cannot be repeated.

Keyword List
------------

**CASE**  *(required)*
A short ASCII label (no spaces) that becomes the name of the FLUNED simulation folder.

**TIME_TREATMENT** *(required)*
Chooses the time integration strategy for the scalar transport solver.

* `steadyState` – runs a pseudo‑steady calculation until convergence.
* `transient`   – using the steady-state velocity field of the input CFD simulation, 
  it calculates the temporal evolution of a the radioisotope concentration field given the 
  chosen boundary conditions
  
**SIMULATION_TYPE** *(required)*
Specify the type of FLUNED simulation.

* `single-isotope` – generates the FLUNED case for a single isotope. Default option if not specified.
* `openmc-multi`   – generates FLUNED cases to study all the relevant radioisotopes. 
  In this mode, an H5 file containing the activation data must be provided. 
  Refer to the user guide documentation for openmc multiple isotopes activation.

**ACTIVATION_FILE** 
In `single-isotope` mode, it is the path to a VTK or VTU file that contains a spatial field of neutron‑activation rate (e.g. N‑16 production). 
Leave the line commented if you prefer a constant or null source term.
In `openmc-multi` mode, it corresponds to the relevant activation data stored in the H5 file generated with the openmc_to_fluned helper function.

**ACTIVATION_DATASET** *(optional, with ACTIVATION_FILE)*
Only for `single-isotope` mode, name of the dataset inside the VTK/VTU file to sample.  Default: `"Value - Total"`.

**ACTIVATION_DATASET_ERROR** *(optional)*
Only for `single-isotope` mode, name of a companion dataset holding MCNP (or other) statistical errors.

**ACTIVATION_CONST**
Only for `single-isotope` mode, uniform volumetric activation rate (SI units 1 s⁻¹ m⁻³). Leave commented if not required.

**ACTIVATION_NORMALIZATION** *(optional)*
Multiplicative factor applied to the sampled activation map before it is written to
the `Source` field.  Keep the default `0` to disable any scaling.

**DECAY_CONSTANT** *(required ≥ 1)*
Radio‑isotope decay constant in s⁻¹. No value is required when running a multi-species simulation using data from OpenMC.

**ISOTOPE** *(optional)*
Isotope string tag, if not specified a `custom` tag is applied and the decay constant must be specified.
No value is required when running a multi-species simulation using data from OpenMC.

**INLET_CONC** *(required)*
Radio‑isotope concentration imposed at all inlet patches \[ m⁻³].
Currently, when running a multi-species simulation using data from OpenMC, the same concentration will applied to all the isotopes. 
This will be improved in future releases.

**MOLECULAR_DIFFUSION** *(required)*
Molecular diffusion coefficient *D* for the species in m² s⁻¹.
Currently, when running a multi-species simulation using data from OpenMC, the same molecular diffusion will applied to all the isotopes. 
This will be improved in future releases.

**SCHMIDT_NUMBER** *(required)* 
Parameter required for the application of the Gradient Diffusion Hypothesis.  
When provided, turbulent diffusivity is taken as *νₜ⁄Sc*.  
A value of `0.7` is often recommended in the literature.
Currently, when running a multi-species simulation using data from OpenMC, the same Schmidt Number will applied to all the isotopes. 
This will be improved in future releases.

**FV_SCHEME** *(optional)*
Controls spatial discretisation accuracy.

* `stable`  (default) – first‑order upwind on convection terms.
* `accurate`          – linearUpwind (blend of upwind and linear).

**CFD_PATH** *(required)*
Absolute or relative path to the **finalised** CFD solution directory.  For OpenFOAM
this is the case root.  For Ansys Fluent supply the folder that holds the
`*.cas.h5`/`*.dat.h5` pair.

**CFD_TYPE** *(required)*
Identifies the CFD format:

* `OpenFoam`        – native OpenFOAM directory tree.
* `fluent-h5-multi` – HDF5 export of a Fluent case in H5 format containing one or more fluid regions. 
  If more than one is present, the `FLUENT_FLUID_REGION_NAME` parameter must be spcified. 
  In this mode, the FLUENT H5 files must contain volume elements with shapes supporrted by OpenFOAM (tetras and hexahedra).

**FLUENT_FLUID_REGION_NAME** 
Exact name of the fluid cell zone to be extracted.  Ignored for OpenFOAM inputs.


