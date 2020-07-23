rule funcotate:
    input: vcf = "{WORKDIR}/VCF/{sample}.passed.vcf",
        ref = config['refGenome']
    output: maf = "{WORKDIR}/VCF/{sample}.funcotator.maf"
    params: funcotatorSources = config["funcotatorSources"],
        genomeVersion = config["genomeVersion"]
    shell:
        "gatk Funcotator " 
        "--variant {input.vcf} " 
        "--reference {input.ref} "
        "--ref-version {params.genomeVersion} "
        "--data-sources-path {params.funcotatorSources} "
        "--output {output.maf} "
        "--output-file-format MAF "
