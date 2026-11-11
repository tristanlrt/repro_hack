rule telecharger_fastq:
    output:
        "data/raw/{srr}.fastq.gz"
    container:
        config["containers"]["sratoolkit"]
    threads: config["threads"]
    shell:
        """
        cd data/raw
        fasterq-dump --threads {threads} --progress {wildcards.srr}
        gzip {wildcards.srr}.fastq
        """


rule telecharger_genome:
    output:
        "data/ref/reference.fasta"
    container:
        config["containers"]["wget"]
    params:
        url=config["reference"]["fasta_url"]
    shell:
        'wget -q -O {output} "{params.url}"'


rule telecharger_annotation:
    output:
        "data/ref/reference.gff"
    container:
        config["containers"]["wget"]
    params:
        url=config["reference"]["gff_url"]
    shell:
        'wget -O {output} "{params.url}"'

rule telecharger_kegg:
    output:
        "data/ref/kegg_{liste}.txt"
    container:
        config["containers"]["wget"]
    params:
        url=lambda wildcards: config["kegg"][wildcards.liste]
    shell:
        'wget -q -O {output} "{params.url}"'
