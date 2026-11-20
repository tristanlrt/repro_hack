# Quantification par transcrit 
# -l A = type de librairie auto (attendu IU) -> lib_format_counts.json
# Taux d'assignation -> aux_info/meta_info.json (percent_mapped)
# Test chr19 = taux ~3-5 %, normal

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
