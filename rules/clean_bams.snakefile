import pandas as pd
from pathlib import Path

REFPATH = config["refGenome"]
REFDICT = config['refGenome'] + ".dict"
WORKDIR = config["outputDir"]
CAPTURE = config["capture"]
BAIT = config["bait"]

include: "bwa_align.snakefile"


def bam_to_clean(wildcards):
    if config["recal_bam"]:
        return (WORKDIR+"/BAM/{sample}.recal.bam", WORKDIR+"/BAM/{sample}.recal.bai")
    return (WORKDIR+"/BAM/{sample}.markdup.bam", WORKDIR+"/BAM/{sample}.markdup.bai")




rule cleaned_bams:
    input:
        bam = bam_to_clean
    output:
        bam = WORKDIR+"/BAM/{sample}.final.bam",
        bai = WORKDIR+"/BAM/{sample}.final.bai"
    shell:
        "mv {input.bam[0]} {output.bam}; mv {input.bam[1]} {output.bai}"



rule markDuplicates:
    input:
        bam = "{WORKDIR}/BAM/{sample}.sorted.bam",
        bai = "{WORKDIR}/BAM/{sample}.sorted.bam.bai"
    output:
        bam = "{WORKDIR}/BAM/{sample}.markdup.bam",
        bai = "{WORKDIR}/BAM/{sample}.markdup.bai",
        metrics = "{WORKDIR}/QC/Picard/{sample}.markdup.metrics"
    params: tmpdir = "{WORKDIR}/BAM",
            ram = config['MarkDuplicates']['ram']
    # shell:	"gatk --java-options '-Xmx{params.ram}G' MarkDuplicates "
    #         "--TMP_DIR={params.tmpdir} "
    #         "--CREATE_INDEX=true -O={output.bam} "
    #         "-I={input.bam} -M={output.metrics}"

    shell: "gatk --java-options '-Xmx{params.ram}G' MarkDuplicates "
            "--TMP_DIR {params.tmpdir} "\
           "--CREATE_INDEX true " \
           "-O {output.bam} -I {input.bam} -M {output.metrics} "

rule genomeDictionary:
	input: REFPATH
	output: REFPATH + ".dict"
	shell: "gatk CreateSequenceDictionary -REFERENCE {input} -OUTPUT {output}"


rule createBaseReacalibrationTable:
    input:
        bai = "{WORKDIR}/BAM/{sample}.markdup.bai",
        ref = REFPATH,
        bam = "{WORKDIR}/BAM/{sample}.markdup.bam",
        dict = REFPATH + ".dict"

    output: "{WORKDIR}/BAM/{sample}.recal.table"
    params: recalVCF1 = config["gnomad"],
        recalVCF2 = config["known_variants"],
        ram = config['MarkDuplicates']['ram']
    shell:
        "gatk IndexFeatureFile -I {params.recalVCF1}; "
        "gatk IndexFeatureFile -I {params.recalVCF2}; "
        "gatk --java-options '-Xmx{params.ram}g' BaseRecalibrator "
        "-R {input.ref} -I {input.bam} --known-sites {params.recalVCF1} "
        "--known-sites {params.recalVCF2} -O {output}"


### Recalibration
rule baseRecalibration:
	input: recaltable = "{WORKDIR}/BAM/{sample}.recal.table",
		ref = REFPATH,
		bam = "{WORKDIR}/BAM/{sample}.markdup.bam",
		bai = "{WORKDIR}/BAM/{sample}.markdup.bai",
		dict = REFDICT
	output:	bam = "{WORKDIR}/BAM/{sample}.recal.bam",
		bai = "{WORKDIR}/BAM/{sample}.recal.bai"
	shell:	"gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr {input.recaltable} "
            "-O {output.bam}"

# QC
rule createPostRecalibrationTable:
    input: ref = REFPATH,
         recalT = "{WORKDIR}/BAM/{sample}.recal.table",
         bam = "{WORKDIR}/BAM/{sample}.final.bam",
         bai = "{WORKDIR}/BAM/{sample}.final.bai",
         dict = REFDICT
    output:	recalTpost = temp("{WORKDIR}/BAM/{sample}.POST.recal.table"),
            plot = "{WORKDIR}/QC/Picard/{sample}.recal_plots.pdf"
    params: recalVCF1 = config["gnomad"],
          recalVCF2 = config["known_variants"],
            ram = config['MarkDuplicates']['ram'],
    shell:	"gatk --java-options '-Xmx{params.ram}g' BaseRecalibrator "
            "-R {input.ref} -I {input.bam} --known-sites {params.recalVCF1} "
              "--known-sites {params.recalVCF2} "
            "-O {output.recalTpost}; "
            "gatk AnalyzeCovariates -before {input.recalT} -after {output.recalTpost} "
            "-plots {output.plot}"

# QC
### PICARD METRICS
rule bedToIntervalFile:
    input: bed = "{bed}",
        dict = REFDICT
    output: "{bed}.IntervalList"
    shell: "gatk BedToIntervalList -INPUT {input.bed} -OUTPUT {output} "
            "--SEQUENCE_DICTIONARY {input.dict}"

rule picardHsMetrics:
    input: bam = "{WORKDIR}/BAM/{sample}.final.bam",
        bai = "{WORKDIR}/BAM/{sample}.final.bai",
        capture = CAPTURE + ".IntervalList",
        ref = REFPATH,
        bait = BAIT + ".IntervalList"
    output:"{WORKDIR}/QC/Picard/{sample}.HsMetrics.txt"
    params: ram=config["MarkDuplicates"]["ram"]
    shell:	"gatk --java-options '-Xmx{params.ram}g' CollectHsMetrics "
            "--INPUT {input.bam} "
            "--OUTPUT {output} "
            "--TARGET_INTERVALS {input.capture} "
            "--BAIT_INTERVALS {input.bait}"

rule picardAlignmentMetrics:
    input: bam = "{WORKDIR}/BAM/{sample}.final.bam",
        bai = "{WORKDIR}/BAM/{sample}.final.bai",
        ref = REFPATH
    output: "{WORKDIR}/QC/Picard/{sample}.AlignmentMetrics.txt"
    shell: "gatk CollectAlignmentSummaryMetrics "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output}"

rule multiqc:
    input:
        expand(WORKDIR+"/QC/Picard/{sample}.HsMetrics.txt", sample=config["samples"].keys()),
        expand(WORKDIR+"/QC/Picard/{sample}.AlignmentMetrics.txt", sample=config["samples"].keys())
    output: "{WORKDIR}/QC/Picard/multiqc_report.html",
        directory("{WORKDIR}/QC/Picard/multiqc_data"),
    shell: "multiqc {WORKDIR}/QC/Picard -o {WORKDIR}/QC/Picard"

# Under Construction
# rule HsMetricsWarnings:
#     input: "{WORKDIR}/QC/Picard/multiqc_data"
#     output: "{WORKDIR}/QC/Picard/WARNINGS.txt"
#     run:
#         warning = ''
#         metrics = config["HsMetrics"]
#         df = pd.read_table(Path(input[0])/"multiqc_picard_HsMetrics.txt", index_col=0)
#         for metric, value in metrics.items():
#             below_metric = ",".join(list(df[df[metric] < value].index))
#             if below_metric:
#                 warning += "{}: {}\n".format(metric, below_metric)
#         with open(output[0], 'w') as fo:
#             fo.write(warning)

