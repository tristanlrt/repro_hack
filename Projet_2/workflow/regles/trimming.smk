# Nettoyage 
# _1P/_2P = paires gardées ; _1U/_2U = orphelines, temporaires
# keepBothReads=True = garde R2 si petit insert (manuel)

rule trimmer:
    input:
        r1="data/raw/{lect}/{srr}_1.fastq.gz",
        r2="data/raw/{lect}/{srr}_2.fastq.gz",
        adaptateurs=config["trimmomatic"]["adaptateurs"]
    output:
        r1p="results/trimmed/{lect}/{srr}_1P.fastq.gz",
        r2p="results/trimmed/{lect}/{srr}_2P.fastq.gz",
        r1u=temp("results/trimmed/{lect}/{srr}_1U.fastq.gz"),
        r2u=temp("results/trimmed/{lect}/{srr}_2U.fastq.gz")
    log:
        "logs/trimmomatic/{lect}/{srr}.log"
    container:
        config["containers"]["trimmomatic"]
    threads: config["threads"]
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
