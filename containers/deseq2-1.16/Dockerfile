FROM mambaorg/micromamba:1.5.8
RUN micromamba install -y -n base -c conda-forge -c bioconda r-base=3.4.1 bioconductor-deseq2=1.16.1 libgfortran=3 && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
