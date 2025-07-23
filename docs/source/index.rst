.. FLUNED documentation master file, created by
   sphinx-quickstart on Wed Jul  9 15:05:01 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


======
FLUNED
======

FLUNED is an open-source, mesh-based tool that reconstructs the time-dependent distribution of activation products in cooling 
circuits exposed to fusion-neutron irradiation.  
Starting from completed CFD run (OpenFOAM or ANSYS Fluent) and the neutron-induced activation data, 
it solves an advection–diffusion–decay problem and outputs concentration fields ready for dose-rate, shielding, or licensing analyses.

Key features
------------

* **CFD-agnostic input** – import meshes and flow fields directly from OpenFOAM/Fluent results; no re-meshing or solver modification required.  
* **Fast parallel solver** – C++ core built on OpenFOAM classes; scales to HPC clusters via native MPI decomposition.  
* **Simple workflow** – `fluned -t` (template), `fluned -i` (pre-process), `FLUNED-solver`, and `fluned-post` streamline the full simulation loop.  
* **Peer-reviewed & validated** – benchmarked against the 2019 FNG water-activation experiment and 
  documented in (`De Pietri et al., 2023 <https://doi.org/10.1016/j.cpc.2023.108807>`_).
* **Quick install** – one-step script (`fluned-repo-install.sh`) for Linux/WSL2, including all Python helper tools.
* **GPL-3.0 License** 


.. toctree::
   :maxdepth: 0
   :caption: Contents:
   :hidden:
   :titlesonly:

  installation
  theory
  usersguide
  validation
  testing
  publications



