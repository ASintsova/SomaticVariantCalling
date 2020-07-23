#configfile: "configs/ancient_config.yaml"
include: "rules/trim.snakefile"
include: "rules/bwa_align.snakefile"
include: "rules/clean_bams.snakefile"
include: "rules/mutect2.smk"
include: "rules/annotate.smk"

rule preprocess:
    input:
        WORKDIR+"/QC/trimmed_reads/multiqc_report.html",
        expand(WORKDIR + "/BAM/{sample}.final.bam", sample=config["samples"].keys()),
        #expand(WORKDIR + "/QC/Picard/{sample}.recal_plots.pdf", sample=config["samples"].keys()),
        WORKDIR + "/QC/Picard/multiqc_report.html"


rule call:
    input:
        #WORKDIR+"/QC/trimmed_reads/multiqc_report.html",
        #WORKDIR + "/QC/Picard/multiqc_report.html",
        expand(WORKDIR +"/VCF/{patientID}.merged.filtered.vcf", patientID=[config["samples"][s]['patientID'] for s in config["samples"].keys()])



rule annotate:
    input:
         expand(WORKDIR + "/VCF/{patientID}.funcotator.maf", patientID=[config["samples"][s]['patientID'] for s in config["samples"].keys()])