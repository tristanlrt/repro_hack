rule deseq2:
    input:
        counts="results/counts/counts.txt",
        samples=config["samples"],
        genes="data/ref/kegg_genes.txt",
        noms="data/ref/aureowiki_noms.tsv",
        ribosome="data/ref/kegg_ribosome.txt",
        trna="data/ref/kegg_trna.txt",
        facteurs="data/ref/kegg_facteurs.txt"
    output:
        table="results/deseq2/resultats_deseq2.tsv",
        resume="results/deseq2/resume.txt",
        ma_complet="results/figures/ma_plot_complet.pdf",
        ma_traduction="results/figures/ma_plot_traduction.pdf"
    container:
        config["containers"]["deseq2"]
    shell:
        "Rscript workflow/scripts/deseq2.R {input.counts} {input.samples} {input.genes} {input.noms} {input.ribosome} {input.trna} {input.facteurs} {output.table} {output.resume} {output.ma_complet} {output.ma_traduction}"
