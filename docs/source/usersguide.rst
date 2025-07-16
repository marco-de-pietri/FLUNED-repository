============
User's Guide
============


This page documents all keywords that can appear in a `FLUNED` input, their
purpose, default behaviour and the range of values they accept.

A generic template can be created by running `fluned -t`

Input Keywords guide
------------------------

Every directive is case‑insensitive and must appear on its own line.  Comment text after a
`#` is ignored.

**CASE**  *(required)*
A short ASCII label (no spaces) that becomes the name of the output folder.  Use it to
distinguish multiple studies that share the same CFD file set.

**`TIME_TREATMENT`** *(required)*
Chooses the time integration strategy for the scalar transport solver.

* `steadyState` – runs a pseudo‑steady calculation until convergence.
* `transient`   – using the steady-state velocity field of the input CFD simulation, 
  it calculates the temporal evolution of a the radioisotope concentration field given the 
  chosen boundary conditions

**`ACTIVATION_FILE`** 
Path to a VTK or VTU file that contains a spatial field of neutron‑activation rate
(e.g. N‑16 production).  Leave the line commented if you prefer a constant or null
source term.

**`ACTIVATION_DATASET`** *(optional, with `ACTIVATION_FILE`)*
Name of the dataset inside the VTK file to sample.  Default: `"Value - Total"`.

**`ACTIVATION_DATASET_ERROR`** *(optional)*
Name of a companion dataset holding MCNP (or other) statistical errors.

**`ACTIVATION_CONST`**
Uniform volumetric activation rate (SI units 1 s⁻¹ m⁻³). Leave commented if not required.

**`ACTIVATION_NORMALIZATION`** *(optional)*
Multiplicative factor applied to the sampled activation map before it is written to
the `Source` field.  Keep the default `0` to disable any scaling.

**`DECAY_CONSTANT`** *(required ≥ 1)*
Radio‑isotope decay constant in **s⁻¹**. No value is required if running a multi-species simulation using data from OpenMC.

**`ISOTOPE`** *(optional)*
String tag (e.g. `N16`) used only for bookkeeping in post‑processing.

**`INLET_CONC`** *(required)*
Radio‑isotope concentration imposed at **all** inlet patches \[ m⁻³].

**`MOLECULAR_DIFFUSION`** *(optional)*
Molecular diffusion coefficient *D* for the species in **m² s⁻¹**.

**`SCHMIDT_NUMBER`**
Dimensionless *Sc = ν⁄D*.  When provided, turbulent diffusivity is taken as
*νₜ⁄Sc*.  Leave the default `0.7` unless you have specific experimental data.

**`FV_SCHEME`** *(optional)*
Controls spatial discretisation accuracy.

* `stable`  (default) – first‑order upwind on convection terms.
* `accurate`          – linearUpwind (blend of upwind and linear).

**`CFD_PATH`** *(required)*
Absolute or relative path to the **finalised** CFD solution directory.  For OpenFOAM
this is the case root.  For Ansys Fluent supply the folder that holds the
`*.cas.h5`/`*.dat.h5` pair.

**`CFD_TYPE`** *(required)*
Identifies the CFD format:

* `OpenFoam`        – native OpenFOAM directory tree.
* `fluent-h5-multi` – HDF5 export of a Fluent case in H5 format containing one or more fluid regions. 
  If more than one is present, the `FLUENT_FLUID_REGION_NAME` parameter must be spcified.

**`FLUENT_FLUID_REGION_NAME`** 
Exact name of the fluid cell zone to be extracted.  Ignored for OpenFOAM inputs.


