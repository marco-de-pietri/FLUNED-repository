#!/bin/bash

# fluned works with python3.11
cd ~
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get -y install python3.11 python3.11-distutils python3.11-venv

curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3.11 get-pip.py
rm get-pip.py

python3.11 -m pip install --user pipx

# Ensure pipx is in PATH for the current script execution
export PATH="$HOME/.local/bin:$PATH"



sudo sh -c "wget -qO - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository -y http://dl.openfoam.org/ubuntu
sudo apt-get update -q
sudo apt-get -y install openfoam12
sudo apt-get -y install libhdf5-dev
sudo apt-get -y install pkg-config

sed -i '$ a . /opt/openfoam12/etc/bashrc' ~/.bashrc

source ~/.bashrc
source /opt/openfoam12/etc/bashrc

cd ~

git clone -v --branch dev --single-branch "https://github.com/marco-de-pietri/FLUNED-repository.git"


cd ~/FLUNED-repository/
pipx install --python python3.11 .


cd ~/FLUNED-repository/FLUNED-solver/
wmake

