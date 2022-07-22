#!/bin/bash

# Download and install `miniconda` in the home directory, initialize and rerun bashrc.
# This script also install mamba for faster dependency solving and installation built in C++ with parallel features.
# It also creates the environment with needed packages for MinPipe.
# Author: Thomaz Guadagnini Ramalheira
# Github: @thzgr
# Email: thzgr@tuta.io

if ! [ -x "$(command -v conda)" ]; then
	HOME_PATH="$(echo $HOME)"
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    bash ~/miniconda.sh -b -p $HOME/miniconda
    
    CONDA_BIN="$HOME_PATH/miniconda/bin/conda"
    
    if [ -f "$CONDA_BIN" ]; then
        cd ~/miniconda/bin
        ./conda init
        source ~/.bashrc
        printf "Conda has been installed and initialized!"
    else
        printf "Conda has not been installed :("
    fi
else
    printf "Conda is already installed!"
fi


if ! [ -x "$(command -v mamba)" ]; then
    conda install mamba -n base -c conda-forge
    printf "Mamba has been installed!"
else
    printf "Mamba is already installed!"

mamba env create -f ../env.yml

mamba init

conda activate minpipe
