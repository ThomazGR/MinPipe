#!/usr/bin/bash

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
conda init
source $HOME/.bashrc

if [ "$CONDA_DEFAULT_ENV" == "" ]; then
    echo "Miniconda has not been installed. Check if there's any errors."
else
    echo "Latest version of Miniconda has been installed, now install the packages!"
fi