# Téléchargement lectures (SRA) et référence (Ensembl 100)

# Lectures paired-end 2x150 pb (LibraryLayout = PAIRED dans le run selector SRA).
# - mode test   : fastq-dump -X <N> ne récupère que les N premières paires
# - mode complet: fasterq-dump (plus rapide) puis compression gzip
# temp() : FASTQ bruts supprimés dès que Trimmomatic les a nettoyés (place disque).
# Pour les garder : ./run.sh --notemp
rule telecharger_fastq:
    output:
        r1=temp("data/raw/{lect}/{srr}_1.fastq.gz"),
        r2=temp("data/raw/{lect}/{srr}_2.fastq.gz")
    container:
        config["containers"]["sratoolkit"]
    threads: config["threads"]
    params:
        max_spots=MAX_SPOTS,
        outdir="data/raw/{lect}"
    shell:
        """
        if [ {params.max_spots} -gt 0 ]; then
            fastq-dump --split-files --gzip -X {params.max_spots} \
                --outdir {params.outdir} {wildcards.srr}
        else
            fasterq-dump --threads {threads} --progress --split-files \
                --outdir {params.outdir} {wildcards.srr}
            gzip {params.outdir}/{wildcards.srr}_1.fastq {params.outdir}/{wildcards.srr}_2.fastq
        fi
        """


# FASTA
rule telecharger_genome:
    priority: 3
    output:
        f"data/ref/{REF}/genome.fa.gz"
    container:
        config["containers"]["wget"]
    params:
        url=GENOME_URL
    shell:
        'wget -q -O {output} "{params.url}"'


# GFF3
rule telecharger_annotation:
    priority: 3
    output:
        f"data/ref/{REF}/annotation.gff3.gz"
    container:
        config["containers"]["wget"]
    params:
        url=GFF_URL
    shell:
        'wget -q -O {output} "{params.url}"'


#  noms des chromosomes
# Conteneur sratoolkit = fournit gzip / grep / sed
rule preparer_reference:
    priority: 3
    input:
        genome=f"data/ref/{REF}/genome.fa.gz",
        gff=f"data/ref/{REF}/annotation.gff3.gz"
    output:
        genome=f"data/ref/{REF}/genome.fa",
        gff=f"data/ref/{REF}/annotation.gff3",
        decoys=f"data/ref/{REF}/decoys.txt"
    container:
        config["containers"]["sratoolkit"]
    shell:
        """
        gunzip -c {input.genome} > {output.genome}
        gunzip -c {input.gff} > {output.gff}
        grep '^>' {output.genome} | cut -d ' ' -f 1 | sed 's/^>//' > {output.decoys}
        """
