#!/bin/bash

# Define variables
INSTALL_DIR="$HOME"

# Navigate to the home directory
cd "$INSTALL_DIR"

# 1) Install required dependencies
echo "Installing required dependencies..."
sudo apt-get update
sudo apt-get -y install software-properties-common git curl wget

# 2) Install Python 3.11 and supporting packages
echo "Adding deadsnakes PPA and installing Python 3.11..."
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get -y install python3.11 python3.11-distutils python3.11-venv

# 3) Install pip for Python 3.11
echo "Installing pip for Python 3.11..."
curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3.11 get-pip.py
rm get-pip.py

# 4) Install pipx using Python 3.11
echo "Installing pipx..."
python3.11 -m pip install --user pipx

# 5) Ensure pipx is in PATH for the script's execution
export PATH="$HOME/.local/bin:$PATH"

# 6) Install OpenFOAM 12
echo "Adding OpenFOAM repository and installing OpenFOAM 12..."
sudo sh -c "wget -qO - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository -y http://dl.openfoam.org/ubuntu
sudo apt-get update -q
sudo apt-get -y install openfoam12

# 7) Install other required dependencies
echo "Installing additional dependencies..."
sudo apt-get -y install libhdf5-dev pkg-config

# 8) Append the OpenFOAM source command to .bashrc
echo "Configuring OpenFOAM environment..."
grep -qxF 'source /opt/openfoam12/etc/bashrc' ~/.bashrc 
source /opt/openfoam12/etc/bashrc  # Source OpenFOAM environment for this subshell


# 9) Clone the FLUNED-repository from GitHub using HTTPS
echo "Cloning FLUNED repository..."
git clone -v --branch master --single-branch "https://github.com/marco-de-pietri/FLUNED-repository.git"

# Navigate to the cloned repository
cd ~/FLUNED-repository/

# 10) Install FLUNED with pipx using Python 3.11
echo "Installing FLUNED using pipx..."
pipx install --python python3.11 .
export PATH="$HOME/.local/bin:$PATH"

# 11) Compile the FLUNED-solver (wmake requires the OpenFOAM environment to be sourced)
echo "Compiling FLUNED-solver..."
cd ~/FLUNED-repository/FLUNED-solver/
wmake

echo "Installation and compilation completed successfully."
