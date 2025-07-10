
============
Installation
============


The recommeded way is to use the installation script present in the ``FLUNED-Repository/installation-script/`` folder.

The tools needs to run in a linux environment, the software requirements are OpenFOAM v12 and the packages listed either in the toml file or in the installation script.

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

