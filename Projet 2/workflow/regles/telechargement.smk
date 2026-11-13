# =============================================================================
# Téléchargement des lectures (SRA) et de la référence (Ensembl release 100)
# =============================================================================

# Lectures paired-end 2x150 pb (LibraryLayout = PAIRED dans le run selector SRA).
# - mode test   : fastq-dump -X <N> ne récupère que les N premières paires
# - mode complet: fasterq-dump (plus rapide) puis compression gzip
rule telecharger_fastq:
    output:
        r1="data/raw/{lect}/{srr}_1.fastq.gz",
        r2="data/raw/{lect}/{srr}_2.fastq.gz"
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


# Génome (FASTA) et annotation (GFF3), compressés
rule telecharger_genome:
    output:
        f"data/ref/{REF}/genome.fa.gz"
    container:
        config["containers"]["wget"]
    params:
        url=GENOME_URL
    shell:
        'wget -q -O {output} "{params.url}"'


rule telecharger_annotation:
    output:
        f"data/ref/{REF}/annotation.gff3.gz"
    container:
        config["containers"]["wget"]
    params:
        url=GFF_URL
    shell:
        'wget -q -O {output} "{params.url}"'


# gffread ne lit pas les fichiers compressés : on décompresse.
# On génère aussi la liste des séquences "leurres" (decoys) pour Salmon :
# ce sont les noms des chromosomes du génome (premier mot des en-têtes FASTA).
# (conteneur sratoolkit utilisé car il fournit gzip/grep/sed, déjà validé au Projet 1)
rule preparer_reference:
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
