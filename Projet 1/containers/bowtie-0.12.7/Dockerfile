FROM mambaorg/micromamba:1.5.8
RUN micromamba install -y -n base -c conda-forge -c bioconda samtools=1.9 wget unzip && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
RUN wget -q -O /tmp/bowtie.zip "https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download" \
 && unzip -q /tmp/bowtie.zip -d /tmp \
 && cp /tmp/bowtie-0.12.7/bowtie /tmp/bowtie-0.12.7/bowtie-build /tmp/bowtie-0.12.7/bowtie-inspect /opt/conda/bin/ \
 && rm -rf /tmp/bowtie.zip /tmp/bowtie-0.12.7
