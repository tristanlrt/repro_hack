#!/bin/bash
if ! command -v conda > /dev/null; then
    export PATH=$HOME/miniforge3/bin:$PATH
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snakemake
snakemake --software-deployment-method apptainer --cores $(nproc) "$@"
