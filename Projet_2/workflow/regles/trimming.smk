# =============================================================================
# Nettoyage des lectures : Trimmomatic v0.39 (mode paired-end)
# =============================================================================
# Sorties : paires conservées (_1P/_2P, utilisées par Salmon) et lectures
# orphelines (_1U/_2U). Toutes temporaires (temp()) : supprimées après la quantification.
# ILLUMINACLIP:<adaptateurs>:2:30:10:2:True = coupe les adaptateurs TruSeq
# (keepBothReads=True conserve la lecture R2 même si elle est le reverse
# complément de R1, recommandé par le manuel pour les petits inserts).

rule trimmer:
    input:
        r1="data/raw/{lect}/{srr}_1.fastq.gz",
        r2="data/raw/{lect}/{srr}_2.fastq.gz",
        adaptateurs=config["trimmomatic"]["adaptateurs"]
    output:
        r1p=temp("results/trimmed/{lect}/{srr}_1P.fastq.gz"),
        r2p=temp("results/trimmed/{lect}/{srr}_2P.fastq.gz"),
        r1u=temp("results/trimmed/{lect}/{srr}_1U.fastq.gz"),
        r2u=temp("results/trimmed/{lect}/{srr}_2U.fastq.gz")
    log:
        "logs/trimmomatic/{lect}/{srr}.log"
    container:
        config["containers"]["trimmomatic"]
    threads: config["threads"]
    priority: 1
    params:
        clip=config["trimmomatic"]["illuminaclip"],
        options=config["trimmomatic"]["options"]
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 \
            {input.r1} {input.r2} \
            {output.r1p} {output.r1u} {output.r2p} {output.r2u} \
            ILLUMINACLIP:{input.adaptateurs}:{params.clip} {params.options} \
            2> {log}
        """
