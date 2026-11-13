# =============================================================================
# Quantification par transcrit : salmon quant (Salmon v1.2.0)
# =============================================================================
# -l A : Salmon détecte le type de librairie (le type retenu est écrit dans
#        lib_format_counts.json ; NEBNext Ultra non directionnel => attendu IU).
# Le taux de lectures assignées est dans aux_info/meta_info.json ("percent_mapped").
# En mode test (chr19 seul) ce taux est faible (~3-5 %) : c'est normal, la
# majorité des lectures viennent des autres chromosomes.

rule quantifier:
    input:
        r1="results/trimmed/{lect}/{srr}_1P.fastq.gz",
        r2="results/trimmed/{lect}/{srr}_2P.fastq.gz",
        index="results/index/salmon_{ref}"
    output:
        quant="results/salmon/{ref}_{lect}/{srr}/quant.sf"
    log:
        "logs/salmon_quant/{ref}_{lect}/{srr}.log"
    container:
        config["containers"]["salmon"]
    threads: config["threads"]
    params:
        libtype=config["salmon"]["libtype"],
        options=config["salmon"]["quant_options"],
        outdir="results/salmon/{ref}_{lect}/{srr}"
    shell:
        """
        salmon quant -i {input.index} -l {params.libtype} \
            -1 {input.r1} -2 {input.r2} \
            {params.options} -p {threads} -o {params.outdir} > {log} 2>&1
        """
