rule indexer_genome:
    input:
        "data/ref/reference.fasta"
    output:
        multiext("results/index/reference", ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt")
    container:
        config["containers"]["bowtie"]
    shell:
        "bowtie-build {input} results/index/reference"
