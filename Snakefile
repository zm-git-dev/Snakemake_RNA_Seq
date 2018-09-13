# Snakemake file for RNA-seq analysis


###############
# Libraries
###############

import os
import pandas as pd
from snakemake.utils import validate, min_version

#############################################
# Configuration and sample sheets
#############################################

configfile: "configs/configs.yaml"
WORKING_DIR         = config["working_dir"]    # where you want to store your intermediate files (this directory will be cleaned up at the end)
RESULT_DIR          = config["result_dir"]      # what you want to keep
# fetch URL to transcriptome multi fasta from configfile
genome_URL = config["refs"]["genome"]
transcriptome_URL = config["refs"]["transcriptome"]

###############
# Helper Functions
###############
def get_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_samples_per_treatment(input_df="units.tsv",colsamples="sample",coltreatment="condition",treatment="control"):
    """This function returns a list of samples that correspond to the same experimental condition"""
    df = pd.read_table(input_df)
    df = df.loc[df[coltreatment] == treatment]
    filtered_samples = df[colsamples].tolist()
    return filtered_samples

##############
# Samples and conditions
##############

units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

SAMPLES = units.index.get_level_values('sample').unique().tolist()

CASES = get_samples_per_treatment(treatment="treatment")
CONTROLS = get_samples_per_treatment(treatment="control")

##############
# Wildcards
##############
wildcard_constraints:
    sample = "[A-Za-z0-9]+"

wildcard_constraints:
    unit = "L[0-9]+"

##############
# Desired output
##############

FASTQC      =   expand(RESULT_DIR + "fastqc/{sample}_{pair}_fastqc.zip", sample=SAMPLES, pair={"forward", "reverse"})
DIFF_EXP    =   RESULT_DIR + "result.csv"

###############
# Final output
################

rule all:
    input:
        FASTQC,
        #DIFF_EXP
    message:
        "Job done!"

###############
# Rules
###############

rule get_genome_fasta:
    output:
        WORKING_DIR + "genome/genome.fasta"
    message:
        "downloading required genomic fasta file"
    shell:
        "wget -O {output} {genome_URL}"

rule get_transcriptome_fasta:
    output:
        WORKING_DIR + "genome/ref_transcriptome.fasta"
    message:
        "downloading required transcriptome fasta file"
    shell:
        "wget -O {output} {transcriptome_URL}"

#rule get_ref_transcriptome_index:
    #input:
    #    WORKING_DIR + "genome/ref_transcriptome.fasta"
    #output:
    #    [WORKING_DIR + "genome/ref_transcriptome.fasta." + i for i in ("nhr", "nin", "nog", "nsd", "nsi", "nsd")]
    #conda:
#        "envs/blast.yaml"
#    shell:
#        "makeblastdb -in {input} -dbtype prot"

rule index:
    input:
        WORKING_DIR + "genome/genome.fasta"
    output:
        [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome with Hisat"

    logs:

    params:
        WORKING_DIR + "genome"
    threads: 10
    shell:
        "hisat2-build --quiet --threads {threads} {input} {params} 2>{log}"

rule trimmomatic:
    input:
        fastq           = get_fastq,
        adapters        = config["adapters"]
    output:
        fw_reads        = WORKING_DIR + "trimmed/{sample}_fw.fq.gz",
        rev_reads       = WORKING_DIR + "trimmed/{sample}_rev.fq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp(WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz")
    message: "trimming {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic/{sample}.log"
    params:
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: 10
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.fastq} "
        "{output.fw_reads} "
        "{output.forwardUnpaired} "
        "{output.rev_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} 2>{log}"

rule fastqc:
    input:
        fwd = WORKING_DIR + "trimmed/{sample}_fw.fq.gz",
        rev = WORKING_DIR + "trimmed/{sample}_rev.fq.gz",
    output:
        fwd = RESULT_DIR + "fastqc/{sample}_forward_fastqc.zip",
        rev = RESULT_DIR + "fastqc/{sample}_reverse_fastqc.zip"
    log:
        "results/logs/fastqc/{sample}.fastqc.log"
    params:
        "results/fastqc/"
    message:
        "Quality check of trimmed samples with FASTQC"
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev} 2>{log}"
#########################################
#########################################
#### FROM HERE I DON'T KNOW KNOW IF IT RUNS
#########################################
#########################################
rule hisat_mapping:
    input:
        fwd = "trimmed/{sample}_fw.fq.gz",
        rev = "trimmed/{sample}_rev.fq.gz",
        ["genome/genome." + str(i) + ".bt2" for i in range(1,5)],
        "genome/genome.rev.1.bt2",
        "genome/genome.rev.2.bt2"
    output:
        bams = expand("mapped/{sample}.bam", sample=samples),
    params:
        index = "genome/genome"
    message:
        "mapping reads to genome to bam files."
    conda:
        "envs/Hisat.yaml"
    shell:
        "python scripts/Hisat.py {params.index_name} {input.fwd} {input.rev} {output.bams}"

rule merge_bams:
    input:
        expand("mapped/{sample}.bam", sample=samples),
    output:
        "mapped/merged.bam"
    conda:
        "envs/Samtools.yaml"
    shell:
        "samtools merge {output} {input}"

rule create_stringtie_transcriptome:
    input:
        "mapped/merged.bam"
    output:
        "genome/stringtie_transcriptome.gtf"
    #params:
        # some parameters
    conda:
        "envs/Stringtie.yaml"
    message:
        "creating transcriptome to stringtie_transcriptome.gtf."
    shell:
        "stringtie {input} {output}"

rule create_counts_table:
    input:
        bams = expand("mapped/{sample}.bam", sample=samples)
        transcriptome = "genome/stringtie_transcriptome.gtf"
    output:
        "results/counts.txt
    #params:
        # some parameters
    conda:
        "subread.yaml"
    shell:
        "python featureCounts.py {input.transcriptome} {output} {inputs.bams}"

rule gtf_to_fasta:
    input:
        tc = "genome/stringtie_transcriptome.gtf"
        genome = "genome/genome.fasta"
    output:
        "genome/stringtie_transcriptome.fasta"
    conda:
        "Tophat.yaml"
    shell:
        "gtf_to_fasta {input.tc} {input.genome} {output}"

rule blast_for_funtions:
    input:
        newTct = "genome/stringtie_transcriptome.fasta"
        refTct = "genome/ref_transcriptome.fasta"
        indexFiles = ["genome/ref_transcriptome.fasta." + i for i in ("nhr", "nin", "nog", "nsd", "nsi", "nsd")]
    output:
        "results/stringtie_transcriptome_blast.txt"
    params:
        evalue = "10"

    conda:
        "envs/blast.yaml"
    shell:
        "blastx -query {input.newTct} -db {input.refTcp} -outfmt 6 qseqid qlen sseqid slen evalue salltitles -out {output} -max_target_seqs 1"

rule DESeq2_analysis:
    input:
        counts    = "results/counts.txt
        functions = "results/stringtie_transcriptome_blast.txt"
    output:
        "results/result.csv"
    message:
        "normalizing read counts en creating differential expression table"
    conda:
        "Deseq.yaml"
    shell:
        "Rscript Deseq.R {input.counts} {input.functions}"
