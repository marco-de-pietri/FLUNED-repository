#!/bin/bash

# Define variables
INSTALL_DIR="$HOME/repos"

# Navigate to the home directory
mkdir -p "$INSTALL_DIR" && cd "$INSTALL_DIR" 

#check for zsh usage
case "$(basename "$SHELL")" in
  zsh)  rc_file="$HOME/.zshrc" ;;
  bash) rc_file="$HOME/.bashrc" ;;
  *)    rc_file="$HOME/.profile" ;;
esac

# 1) Install required dependencies
echo "Installing required dependencies..."
sudo apt-get update
sudo apt-get -y install software-properties-common git curl wget


# 2) Install OpenFOAM 12
echo "Adding OpenFOAM repository and installing OpenFOAM 12..."
sudo sh -c "wget -qO - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository -y http://dl.openfoam.org/ubuntu
sudo apt-get update -q
sudo apt-get -y install openfoam12

# 3) Install other required dependencies
echo "Installing additional dependencies..."
sudo apt-get -y install libhdf5-dev pkg-config

# 4) Append the OpenFOAM source command to .bashrc
echo "Configuring OpenFOAM environment..."
grep -qxF 'source /opt/openfoam12/etc/bashrc' "$rc_file" \
  || printf '\nsource /opt/openfoam12/etc/bashrc\n' >> "$rc_file"
source /opt/openfoam12/etc/bashrc  # Source OpenFOAM environment for this subshell


# 5) Clone the FLUNED-repository from GitHub using HTTPS
echo "Cloning FLUNED repository..."
DIR="FLUNED-repository"
if [ -d "$DIR" ]; then
  rm -rf "$DIR"
fi
git clone -v --branch master --single-branch "https://github.com/marco-de-pietri/FLUNED-repository.git" "$DIR"

# Navigate to the cloned repository
cd "$INSTALL_DIR/FLUNED-repository"

#  Install uv
echo "Installing required dependencies..."
wget -qO- https://astral.sh/uv/install.sh | sh

#  Ensure uv is in PATH for the script's execution
export PATH="$HOME/.local/bin:$PATH"

uv tool install . --python 3.13



# 6) Install FLUNED with pipx using Python 3.11
echo "Installing FLUNED using pipx..."
pipx install --python python3.11 .

grep -qxF 'export PATH="$HOME/.local/bin:$PATH"' "$rc_file" \
  || printf '\nexport PATH="$HOME/.local/bin:$PATH"\n' >> "$rc_file"



# 7) Compile the FLUNED-solver (wmake requires the OpenFOAM environment to be sourced)
echo "Compiling FLUNED-solver..."
cd "$INSTALL_DIR/FLUNED-repository/FLUNED-solver"
wmake

echo "Installation and compilation completed successfully."

