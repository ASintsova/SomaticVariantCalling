WORKDIR = config["outputDir"]

def getFastPair(wildcards):
    forward = config["samples"][wildcards.sample]["fwd"]
    reverse = config["samples"][wildcards.sample]["rvr"]
    return forward, reverse

rule trimmed:
    input:  WORKDIR+"/QC/trimmed_reads/multiqc_report.html"


rule trim:
    input:
        pair = getFastPair
    params:
        q = config["skewer"]["q"],
        Q = config["skewer"]["Q"],
        l = config["skewer"]["l"],
        x = config["skewer"]["x"],
        y = config["skewer"]["y"],
        prefix = "{WORKDIR}/trimmed_reads/{sample}",
    log:
        "{WORKDIR}/logs/trimmed_reads/{sample}.log"
    output:
        temp("{WORKDIR}/trimmed_reads/{sample}-trimmed-pair1.fastq"),
        temp("{WORKDIR}/trimmed_reads/{sample}-trimmed-pair2.fastq")
    shell:
        "mkdir -p {WORKDIR}/trimmed_reads; "
        "skewer -x {params.x} "
        "-y {params.y} "
        "-Q {params.Q} "
        "-q {params.q} "
        "-l {params.l} "
        "{input} "
        "-o {params.prefix} 1> {log} 2> {log}"


rule fastqc:
    input: "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair1.fastq",
        "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair2.fastq"
    params: outdir = WORKDIR+"/QC/trimmed_reads"
    output: "{WORKDIR}/QC/trimmed_reads/{sample}-trimmed-pair1_fastqc.html",
        "{WORKDIR}/QC/trimmed_reads/{sample}-trimmed-pair2_fastqc.html"
    shell: "fastqc -o {params.outdir} {input}; "


rule multiqcOnFastqc:
    input:
        expand(WORKDIR+"/QC/trimmed_reads/{sample}-trimmed-pair1_fastqc.html", sample=config["samples"].keys()),
        expand(WORKDIR+"/QC/trimmed_reads/{sample}-trimmed-pair2_fastqc.html", sample=config["samples"].keys())
    output:
        WORKDIR+"/QC/trimmed_reads/multiqc_report.html"
    params: outdir=WORKDIR+"/QC/trimmed_reads"
    shell: "multiqc {params.outdir} -o {params.outdir}"

