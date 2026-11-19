rule mapper:
    input:
        fastq="results/trimmed/{srr}.trimmed.fastq.gz",
        index=multiext("results/index/reference", ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt")
    output:
        bam="results/bam/{srr}.bam",
        bai="results/bam/{srr}.bam.bai"
    container:
        config["containers"]["bowtie"]
    threads: config["threads"]
    shell:
        """
        bowtie -p {threads} -S results/index/reference <(gunzip -c {input.fastq}) | samtools sort -@ {threads} > {output.bam}
        samtools index {output.bam}
        """
