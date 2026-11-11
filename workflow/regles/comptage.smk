rule compter:
    input:
        bam=expand("results/bam/{srr}.bam", srr=SAMPLES),
        bai=expand("results/bam/{srr}.bam.bai", srr=SAMPLES),
        gff="data/ref/reference.gff"
    output:
        "results/counts/counts.txt"
    container:
        config["containers"]["subread"]
    threads: config["threads"]
    shell:
        "featureCounts -t gene -g ID -s 1 -F GTF -T {threads} -a {input.gff} -o {output} {input.bam}"
