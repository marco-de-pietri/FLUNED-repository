#!/bin/bash

# Ensure the script is being sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Please run this script using 'source install_fluned.sh' or '. install_fluned.sh'"
    exit 1
fi

# Exit immediately if a command exits with a non-zero status
set -e

# Define variables
INSTALL_DIR="$HOME"

# Navigate to the home directory
cd "$INSTALL_DIR"

# 1) Install required dependencies
sudo apt-get update
sudo apt-get -y install software-properties-common git curl wget

# 2) Install Python 3.11 and supporting packages
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get -y install python3.11 python3.11-distutils python3.11-venv

# 3) Install pip for Python 3.11
curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3.11 get-pip.py
rm get-pip.py

# 4) Install pipx using Python 3.11
python3.11 -m pip install --user pipx

# 5) Ensure pipx is in PATH for the current shell
export PATH="$HOME/.local/bin:$PATH"

# 6) Install OpenFOAM 12
sudo sh -c "wget -qO - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository -y http://dl.openfoam.org/ubuntu
sudo apt-get update -q
sudo apt-get -y install openfoam12

# 7) Install other required dependencies
sudo apt-get -y install libhdf5-dev pkg-config

# 8) Source OpenFOAM environment directly within the script
source /opt/openfoam12/etc/bashrc

# 9) (Optional) Append the OpenFOAM source command to .bashrc for future sessions
grep -qxF 'source /opt/openfoam12/etc/bashrc' ~/.bashrc || echo 'source /opt/openfoam12/etc/bashrc' >> ~/.bashrc

# 10) Clone the FLUNED-repository from GitHub using HTTPS
git clone -v --branch dev --single-branch "https://github.com/marco-de-pietri/FLUNED-repository.git"

# Navigate to the cloned repository
cd ~/FLUNED-repository/

# 11) Install FLUNED with pipx using Python 3.11
pipx install --python python3.11 .

# 12) Compile the FLUNED-solver (wmake requires the OpenFOAM environment to be sourced)
cd ~/FLUNED-repository/FLUNED-solver/
wmake

# 13) Update PATH for pipx-installed modules immediately
export PATH="$HOME/.local/bin:$PATH"


