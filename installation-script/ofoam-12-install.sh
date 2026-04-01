#check for zsh usage
case "$(basename "$SHELL")" in
zsh) rc_file="$HOME/.zshrc" ;;
bash) rc_file="$HOME/.bashrc" ;;
*) rc_file="$HOME/.profile" ;;
esac

# Install required dependencies
echo "Installing required dependencies..."
sudo apt-get update
sudo apt-get -y install software-properties-common git curl wget

# Install OpenFOAM 12
echo "Adding OpenFOAM repository and installing OpenFOAM 12..."
sudo sh -c "wget -qO - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository -y http://dl.openfoam.org/ubuntu
sudo apt-get update -q
sudo apt-get -y install openfoam12

# Install other required dependencies
echo "Installing additional dependencies..."
sudo apt-get -y install libhdf5-dev pkg-config

# Append the OpenFOAM source command to the shell rc file
echo "Configuring OpenFOAM environment..."
foam_line='[ -r /opt/openfoam12/etc/bashrc ] && . /opt/openfoam12/etc/bashrc'
grep -qxF "$foam_line" "$rc_file" ||
  printf '\n# OpenFOAM\n%s\n' "$foam_line" >>"$rc_file"

source /opt/openfoam12/etc/bashrc
