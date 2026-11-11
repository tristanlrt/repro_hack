FROM mambaorg/micromamba:1.5.8
RUN micromamba install -y -n base -c conda-forge -c bioconda cutadapt=1.11 && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
