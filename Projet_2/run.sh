#!/bin/bash
if ! command -v conda > /dev/null; then
    export PATH=$HOME/miniforge3/bin:$PATH
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snakemake
# Images hors du projet : Snakemake ne met pas leur chemin entre guillemets, donc un espace (ex. un dossier parent avec espace) casse apptainer
snakemake --software-deployment-method apptainer --apptainer-prefix "$HOME/.cache/apptainer-images" --cores $(nproc) "$@"
