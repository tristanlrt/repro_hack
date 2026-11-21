rule trimmer:
    input:
        "data/raw/{srr}.fastq.gz"
    output:
        "results/trimmed/{srr}.trimmed.fastq.gz"
    container:
        config["containers"]["cutadapt"]
    shell:
        "cutadapt -q 20 -m 25 -o {output} {input}"
