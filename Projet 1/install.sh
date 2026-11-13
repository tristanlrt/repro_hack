#!/bin/bash
set -e

if ! command -v apptainer > /dev/null; then
    wget -q -O /tmp/apptainer.deb https://github.com/apptainer/apptainer/releases/download/v1.5.5/apptainer_1.5.5_amd64.deb
    sudo apt-get update
    sudo apt-get install -y /tmp/apptainer.deb
    rm /tmp/apptainer.deb
fi

if ! command -v conda > /dev/null; then
    wget -q -O /tmp/miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash /tmp/miniforge.sh -b -p $HOME/miniforge3
    rm /tmp/miniforge.sh
    export PATH=$HOME/miniforge3/bin:$PATH
fi

source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | grep -q "^snakemake "; then
    conda create -y -n snakemake -c conda-forge -c bioconda snakemake=9.27.1 pandas
fi

apptainer --version
conda activate snakemake
snakemake --version
