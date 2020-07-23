
REFPATH = config["refGenome"]
WORKDIR = config["outputDir"]

include: "trim.snakefile"

def get_sample_pateint(wildcards):
    return config["samples"][wildcards.sample]["patientID"]

def get_sample_tissue(wildcards):
    return config["samples"][wildcards.sample]["tissue"]

rule bwaIndex:
    input: REFPATH
    output: REFPATH + ".bwt"
    shell: "bwa index {input}"

rule aligned:
    input:
        expand(WORKDIR+"/BAM/{sample}.sorted.bam", sample=config["samples"].keys())

rule align:
    input:
        ref = REFPATH,
        frwd = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair1.fastq",
        rvr = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair2.fastq",
        index = REFPATH + ".bwt"
    output:	bam = temp("{WORKDIR}/BAM/{sample}.sorted.bam"),
        #bai = temp("{WORKDIR}/BAM/{sample}.sorted.bam.bai")
        bai = "{WORKDIR}/BAM/{sample}.sorted.bam.bai"
    params: pid = get_sample_pateint,
        tissue = get_sample_tissue,
        #threads = config['threads']
    shell:	"bwa mem  -t 4 "
            "-R '@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tSM:{params.pid}_{params.tissue}\\tPL:Illumina\\tCN:CENTER' "
            "-M {input.ref} "
            "{input.frwd} {input.rvr} | samtools sort --reference {input.ref} "
            "-l 9 -@ 4 -O bam -o {output.bam} -T {output.bam}.tmp; samtools index {output.bam}"