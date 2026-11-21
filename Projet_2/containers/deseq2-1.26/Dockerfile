# DESeq2 1.26.0 (Bioconductor 3.10) sous R 3.6.1
# Recette équivalente à l'image Biocontainers utilisée dans config.yaml
# (bioconductor-deseq2:1.26.0--r36he1b5a44_0), même principe que le Projet 1.
FROM mambaorg/micromamba:1.5.8
RUN micromamba install -y -n base -c conda-forge -c bioconda r-base=3.6.1 bioconductor-deseq2=1.26.0 && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
