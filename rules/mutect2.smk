import pandas as pd
#configfile: "configs/test_config.yaml"
WORKDIR = config['outputDir']

#samples = pd.read_csv("configs/test_sample.csv")


# Get optional params for Mutect2 out of the config

def optional_params():
    extra = ''
    extra += '--germline-resource {} '.format(config['gnomad'])if config['gnomad'] else ''
    extra += "-pon {} ".format(config['pon']) if config['pon'] else ''
    for param in config['parameters']:
        if config['parameters'][param]:
            extra += '-{} {} '.format(param, config['parameters'][param])
    return extra


def patient2tumor(wildcards):
    if config['preprocess']:

        for sample, values in config['samples'].items():
            if (str(values['patientID']) == wildcards.patientID) and (values['tissue'] == 'tumor'):
                return "{}/BAM/{}.final.bam".format(WORKDIR, sample)
        return "No tumor sample found in config"
    else:
        for sample, values in config['samples'].items():
            if (str(values['patientID']) == wildcards.patientID) and (values['tissue'] == 'tumor'):
                return values['bam']


def allTumorBams(wildcards):
    samples = []
    if config['preprocess']:
        for sample, values in config['samples'].items():
            if values['tissue'] == 'tumor':
                samples.append(sample)
        return expand(WORKDIR + "/BAM/{sample}.final.bam", sample = samples)
    else:
        for sample, values in config['samples'].items():
            if values['tissue'] == 'tumor':
                samples.append(values['bam'])
        return samples
    #return "No tumor sample found in config"


def patient2normal(wildcards):
    if config['preprocess']:
        for sample, values in config['samples'].items():
            if (str(values['patientID']) == wildcards.patientID) and (values['tissue'] == 'normal'):
                return "{}/BAM/{}.final.bam".format(WORKDIR, sample)
        return "No tumor sample found in config"
    else:
        for sample, values in config['samples'].items():
            if (str(values['patientID']) == wildcards.patientID) and (values['tissue'] == 'normal'):
                return values['bam']



rule all:
    input: WORKDIR + "/VCF/p1.pileupsummaries.table"


# Split Intervals
rule Mutect2IntervalSplit:
    input:
        capture_intervals = config['capture'],
        ref = config['refGenome']
    output: temp(expand(WORKDIR+"/scatter/{sample}-scattered.interval_list", sample=[str(i).zfill(4) for i in range(0, config['scatter'])]))
    params: scatter = config['scatter']
    shell: "gatk SplitIntervals "
           "-R {input.ref} "
           "-L {input.capture_intervals} "
           "--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION "
           "--scatter-count {params.scatter} "
           "-O {WORKDIR}/scatter" # todo make this a tempdir



# Variant Calling
rule mutect2:
    input:
        tumor = patient2tumor,
        normal = patient2normal,
        ref = config['refGenome'],
        reffai = config['refGenome'] + ".fai",
        refdict = config['refGenome'].rstrip(".fasta") + ".dict",
        interval = "{WORKDIR}/scatter/{scatter}-scattered.interval_list"
    output:
        vcf = temp("{WORKDIR}/VCF/{patientID}-{scatter}.vcf"),
        stat = temp("{WORKDIR}/VCF/{patientID}-{scatter}.vcf.stats"),
        i = temp("{WORKDIR}/VCF/{patientID}-{scatter}.vcf.idx")
    params:
        extra = optional_params()
    shell: "mkdir -p {WORKDIR}/VCF; "
        "gatk Mutect2 "
        "-R {input.ref} "
        "-I {input.tumor} "
        "-I {input.normal} "
        "{params.extra} "
        "-L {input.interval} "
        "-normal $(samtools view -H {input.normal} |grep 'RG'|head -1|awk -F 'SM:' '{{print $2}}'|cut -f1) "
        "-O {output.vcf}"


# Merge Interval VCFs

rule MergeVCFs:
    input:
        vcfs = list(map(lambda x: "{WORKDIR}/VCF/{patientID}" + x, ["-{}.vcf".format(s) for s in [str(i).zfill(4) for i in range(0, config['scatter'])]])),
        stats = list(map(lambda x: "{WORKDIR}/VCF/{patientID}" + x, ["-{}.vcf.stats".format(s) for s in [str(i).zfill(4) for i in range(0, config['scatter'])]])),
        indexes = list(map(lambda x: "{WORKDIR}/VCF/{patientID}" + x, ["-{}.vcf.idx".format(s) for s in [str(i).zfill(4) for i in range(0, config['scatter'])]])),
    output: mvcf = "{WORKDIR}/VCF/{patientID}.unfiltered.vcf",
           mstats =  "{WORKDIR}/VCF/{patientID}.unfiltered.vcf.stats"
    params: vcf_files = lambda wildcards, input: ' '.join(["-I {} ".format(f) for f in input.vcfs]),
            stat_files = lambda wildcards, input: ' '.join(["-stats {}.stats ".format(f) for f in input.vcfs])
    shell: "gatk MergeVcfs {params.vcf_files} -O {output.mvcf}; "
           "gatk MergeMutectStats {params.stat_files} -O {output.mstats}"




# Potentially has a bug,  check for a version of GATK where the bug is fixed...


# Build Read orientation model

rule mutect2F1R2:
    input:
        tumor_bams = allTumorBams,
        ref = config['refGenome'],
        reffai =  config['refGenome'] + ".fai",
        interval = "{WORKDIR}/scatter/{scatter}-scattered.interval_list"
    params:
        extra = optional_params(),
        files = lambda wildcards, input: ' '.join(["-I {}".format(f) for f in input.tumor_bams])
    output: temp("{WORKDIR}/VCF/f1r2-{scatter}.tar.gz"),
        temp("{WORKDIR}/VCF/{scatter}-F1R2.unfiltered.vcf"),
        temp("{WORKDIR}/VCF/{scatter}-F1R2.unfiltered.vcf.idx"),
        temp("{WORKDIR}/VCF/{scatter}-F1R2.unfiltered.vcf.stats")
    shell: "mkdir -p {WORKDIR}/VCF; "
        "gatk Mutect2 -R {input.ref} "
        "{params.files} "
        "{params.extra} "
        "-L {input.interval} "
        "--f1r2-tar-gz {output[0]} -O {output[1]}"

rule learnReadOrientationModel:
    input:
        f1r2 = expand(WORKDIR+"/VCF/f1r2-{scatter}.tar.gz", scatter=[str(i).zfill(4) for i in range(config['scatter'])])
    output:
        "{WORKDIR}/VCF/read-orientation-model.tar.gz"
    params:
        f1r2 = lambda wildcards, input: ' '.join([" -I {} ".format(f) for f in input.f1r2])
    shell: "gatk LearnReadOrientationModel {params.f1r2} -O {output}"


# Estimate Contamination

rule GetPileupSummaries:
    input:
        tumor = patient2tumor
    output:
        psum = "{WORKDIR}/VCF/{patientID}.pileupsummaries.table"
    params:
        variants = config['known_variants']
    shell:
        "gatk GetPileupSummaries -I {input.tumor} "
        "-V {params.variants} "
        "-L {params.variants} "
        "-O {output.psum}"


rule CalculateContamination:
    input:
        psum =  "{WORKDIR}/VCF/{patientID}.pileupsummaries.table"
    output:
        segtable = "{WORKDIR}/VCF/{patientID}.segmentation.table",
        conttable = "{WORKDIR}/VCF/{patientID}.contamination.table"
    shell:
        "gatk CalculateContamination "
        "-I {input.psum} "
        "-tumor-segmentation {output.segtable} "
        "-O {output.conttable} "




# Filter Variants

rule filterVariants:
    input:
        ref = config['refGenome'],
        vcf = "{WORKDIR}/VCF/{patientID}.unfiltered.vcf",
        ROmodel = "{WORKDIR}/VCF/read-orientation-model.tar.gz",
        segtable = "{WORKDIR}/VCF/{patientID}.segmentation.table",
        conttable = "{WORKDIR}/VCF/{patientID}.contamination.table"
    output: "{WORKDIR}/VCF/{patientID}.merged.filtered.vcf"
    shell: "gatk FilterMutectCalls "
           "-R {input.ref} "
		   "-V {input.vcf} "
           "--tumor-segmentation {input.segtable} "
           "--contamination-table {input.conttable} "
           "--ob-priors {input.ROmodel} "
           "-O {output}"


rule selectPassed:
    input: vcf = "{WORKDIR}/VCF/{sample}.merged.filtered.vcf"
    output: "{WORKDIR}/VCF/{sample}.passed.vcf"
    shell: "gatk SelectVariants -V {input.vcf} --exclude-filtered  -O {output}"



# Helper Rules


rule fasta_index:
    input: "{sample}.fasta"
    output: "{sample}.fasta.fai"
    shell: "samtools faidx {input}"


rule fasta_dict:
    input: "{sample}.fasta"
    output:"{sample}.dict"
    shell: "gatk CreateSequenceDictionary -R {input}"












#
#
# #### Create panel of normals
#
# rule mutect2forPON:
#     input: bams = patient2normal,
#         ref = REFPATH,
#         reffai = REFFAI,
#         interval = "{WORKDIR}/scatter/{scatter}-scattered.interval_list"
#     params:
#         germRes = germRes,
#     output: temp("{WORKDIR}/PON/PON-{patientID}-{scatter}.vcf")
#     shell: "mkdir -p {WORKDIR}/PON; gatk Mutect2 "
# 		   "-R {input.ref} "
# 		   "-I {input.bams} -L {input.interval} {params.germRes} "
# 		   "--max-mnp-distance 0 "
# 		   "-O {output}"
#
#
#
# # rule mutect2forPON:
# #     input: bams = patient2normal,
# #         ref = REFPATH,
# #         reffai = REFFAI
# #     params:
# #         intList = intList,
# #         germRes = germRes,
# #     output: "{WORKDIR}/PON/PON-{patientID}.vcf"
# #     shell: "mkdir -p {WORKDIR}/PON; gatk Mutect2 "
# # 		   "-R {input.ref} "
# # 		   "-I {input.bams} {params.intList} {params.germRes} "
# # 		   "--max-mnp-distance 0 "
# # 		   "-O {output}"
#
# rule callPON:
#     input: normal_vcfs = expand(WORKDIR+"/PON/PON-{samples}.merged.vcf", samples=config['patients'].keys()),
#         ref = REFPATH,
#         interval = "{WORKDIR}/temp/{scatter}-scattered.interval_list"
#     output: vcf = "{WORKDIR}/PON/PON-{scatter}.vcf"
#     params: files = lambda wildcards, input: ' '.join(["-V {}".format(f) for f in input.normal_vcfs]),
#         germRes = germRes,
#         pondb = "{WORKDIR}/PON/pon_db "
#     shell: "gatk GenomicsDBImport "
#        "-R {input.ref} "
#        "-L {input.interval} "
#        "--genomicsdb-workspace-path {params.pondb} {params.files}; "
#        "gatk CreateSomaticPanelOfNormals "
#        "-R {input.ref} "
#        "{params.germRes} "
#        "-V gendb://{params.pondb} "
#        "-O {output.vcf}; rm -r {params.pondb} "
#

#
#