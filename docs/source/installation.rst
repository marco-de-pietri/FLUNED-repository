============
Installation
============


The recommeded installation method is running the  installation script present in the ``FLUNED-Repository/installation-script/`` folder.

This scripts installs the required dependencies, compiles the FLUNED solver and install the python modules which expose the pre- and post-processor scripts.

The tools needs to run in a linux environment (linux of WSL2 in windows). FLUNED requires the installation of OpenFOAM v12, on which it is based. 
The installation script uses the uv python package manager for the python supporting scripts.

The simplest way is to download the script and run it with the following commands:

1. Download

   .. code-block:: bash

      wget https://raw.githubusercontent.com/marco-de-pietri/FLUNED-repository/refs/heads/master/installation-script/fluned-repo-install.sh

2. Make it executable

   .. code-block:: bash

      chmod +x fluned-repo-install.sh

3. Install it

   .. code-block:: bash

      ./fluned-repo-install.sh

4. Restart shell or source the .bashrc file

   .. code-block:: bash

      source ~/.bashrc

Users can also check the script and, if required, install it with a different workflow or tooling.

