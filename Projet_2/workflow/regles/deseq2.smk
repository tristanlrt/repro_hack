# =============================================================================
# Analyse différentielle par transcrit : DESeq2 1.26.0 (R 3.6)
# =============================================================================
# L'article écrit « DEseq2 package (v.3.6.1) », qui n'est pas une version de
# DESeq2. Hypothèse : c'est la version de R (3.6.1). Bioconductor 3.10 est la
# version sortie pour R 3.6.1 ; elle contient DESeq2 1.26.0.
# Conteneur : image Biocontainers bioconductor-deseq2 1.26.0 (référencée par
# Damona). Recette équivalente : containers/deseq2-1.26/Dockerfile.
#
# Entrée : les 12 quant.sf d'un même mode (test ou complet).
# Sorties (dans results/deseq2/<ref>_<lect>/) :
#   comptes_transcrits.tsv   matrice de comptes (NumReads arrondis)
#   resultats_d11.tsv / _d33 résultats DESeq2 complets + colonne DET (up/down/non)
#   comparaison_article.tsv  nombre de DETs pour 4 variantes vs l'article
#   venn_dets.pdf            diagramme de Venn (Fig. 2B)

rule deseq2:
    input:
        samples=config["samples"],
        script="workflow/scripts/deseq2.R",
        quant=expand("results/salmon/{{ref}}_{{lect}}/{srr}/quant.sf", srr=SAMPLES)
    output:
        comptes="results/deseq2/{ref}_{lect}/comptes_transcrits.tsv",
        d11="results/deseq2/{ref}_{lect}/resultats_d11.tsv",
        d33="results/deseq2/{ref}_{lect}/resultats_d33.tsv",
        comparaison="results/deseq2/{ref}_{lect}/comparaison_article.tsv",
        venn="results/deseq2/{ref}_{lect}/venn_dets.pdf"
    log:
        "logs/deseq2/{ref}_{lect}.log"
    container:
        config["containers"]["deseq2"]
    params:
        modele=config["deseq2"]["modele"],
        colonne_p=config["deseq2"]["colonne_p"],
        seuil_p=config["deseq2"]["seuil_p"],
        seuil_lfc=config["deseq2"]["seuil_log2fc"]
    shell:
        """
        Rscript {input.script} {input.samples} \
            {output.comptes} {output.d11} {output.d33} {output.comparaison} {output.venn} \
            {params.modele} {params.colonne_p} {params.seuil_p} {params.seuil_lfc} \
            {input.quant} > {log} 2>&1
        """
