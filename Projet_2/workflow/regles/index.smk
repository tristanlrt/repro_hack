# Transcriptome (gffread) + index Salmon 1.2.0

# Transcrits = exons concaténés (gffread -w), identifiants ENSMUST
rule construire_transcriptome:
    input:
        genome=f"data/ref/{REF}/genome.fa",
        gff=f"data/ref/{REF}/annotation.gff3"
    output:
        f"data/ref/{REF}/transcriptome.fa"
    container:
        config["containers"]["gffread"]
    shell:
        "gffread -w {output} -g {input.genome} {input.gff}"


# Index avec decoys : transcrits + génome (gentrome), génome = leurre
# decoys: false -> transcrits seuls
rule indexer_salmon:
    input:
        transcrits=f"data/ref/{REF}/transcriptome.fa",
        genome=f"data/ref/{REF}/genome.fa",
        decoys=f"data/ref/{REF}/decoys.txt"
    output:
        directory(f"results/index/salmon_{REF}")
    log:
        f"logs/salmon_index/{REF}.log"
    container:
        config["containers"]["salmon"]
    threads: config["threads"]
    params:
        kmer=config["salmon"]["kmer"],
        decoys="true" if config["salmon"]["decoys"] else "false",
        gentrome=f"data/ref/{REF}/gentrome.fa"
    shell:
        """
        if [ {params.decoys} = "true" ]; then
            cat {input.transcrits} {input.genome} > {params.gentrome}
            salmon index -t {params.gentrome} -d {input.decoys} \
                -i {output} -k {params.kmer} -p {threads} > {log} 2>&1
            rm {params.gentrome}
        else
            salmon index -t {input.transcrits} \
                -i {output} -k {params.kmer} -p {threads} > {log} 2>&1
        fi
        """
