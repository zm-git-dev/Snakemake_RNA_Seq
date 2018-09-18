configfile: "configs/configs.yaml"

# fetch URL to transcriptome multi fasta from configfile
genome_URL = config["refs"]["genome"]
#transcriptome_fasta_URL = config["refs"]["transcriptomeFas"]
transcriptome_gtf_URL = config["refs"]["transcriptomeGtf"]
# create lists containing samplenames, conditions, path/filenames of the fw-reads(R1)
# and the path/filenames of the rev reads(R2) from the file: data/sampls.txt
import pandas as pd
SAMPLES = list(pd.read_table(config["units"])["sample"])
#conditions = list(pd.read_table("data/samples.txt")["condition"])
#R1 = list(pd.read_table("data/samples.txt")["fq1"])
#R2 = list(pd.read_table("data/samples.txt")["fq2"])
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)
#SAMPLES = units.index.get_level_values('sample').unique().tolist()

rule all:
    input:
        countsTable = "results/counts.txt",
        #resultsFile ="results/result.csv",
        fwd = expand("results/fastqc/{sample}_fw_fastqc.zip", sample = SAMPLES),
        rev = expand("results/fastqc/{sample}_rev_fastqc.zip", sample = SAMPLES)
    message:
        "Job done!"

rule get_genome_fasta:
    output:
        "genome/genome.fasta"
    message:"downloading required genomic fasta file"
    shell: "wget -O {output} {genome_URL}"

rule get_transcriptome_fasta:
    output:
        "genome/ref_transcriptome.fasta"
    message:"downloading required transcriptome fasta file"
    shell: "wget -O {output} {transcriptome_fasta_URL}"

rule get_transcriptome_gtf:
    output:
        "genome/ref_transcriptome.gff"
    message:"downloading required transcriptome gtf file"
    shell: "wget -O {output} {transcriptome_gtf_URL}"

rule get_ref_transcriptome_index:
    input:
        "genome/ref_transcriptome.fasta"
    output:
        ["genome/ref_transcriptome.fasta." + i for i in ("nhr", "nin", "nog", "nsd", "nsi")]
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

rule trimmomatic:
    input:
        fq1 = "data/{SAMPLES}_R1.sub.fq",
        fq2 = "data/{SAMPLES}_R2.sub.fq",
        adapters = config["adapters"]
    output:
        fw_reads = "trimmed/{SAMPLES}_fw.fq",
        rev_reads = "trimmed/{SAMPLES}_rev.fq",
        forwardUnpaired = "trimmed/{SAMPLES}_forward_unpaired.fastq",
        reverseUnpaired = "trimmed/{SAMPLES}_reverse_unpaired.fastq"
#    message: "trimming reads"
#        "logs/trimmomatic/{SAMPLES}.log"
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
    threads: 1
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.fq1} "
        "{input.fq2} "
        "{output.fw_reads} "
        "{output.forwardUnpaired} "
        "{output.rev_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen}" #" 2>{log}"

rule fastqc:
    input:
        fwd = "trimmed/{SAMPLES}_fw.fq",
        rev = "trimmed/{SAMPLES}_rev.fq",
    output:
        fwd="results/fastqc/{SAMPLES}_fw_fastqc.zip",
        rev="results/fastqc/{SAMPLES}_rev_fastqc.zip"
    log:
        "results/logs/fastqc/{SAMPLES}.fastqc.log"
    params:
        "results/fastqc/"
    message:
        "Quality check of trimmed samples with FASTQC" 		#removed, it was not working
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev}"

rule index:
    input:
        "genome/genome.fasta"
    output:
        ["genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    message:"indexing genome"
    params:
        "genome/genome"
    conda:
        "envs/Hisat.yaml"
    threads: 10
    shell:"hisat2-build -p {threads} {input} {params}"

rule hisat_mapping:
    input:
        fwd   = "trimmed/{SAMPLES}_fw.fq",
        rev   = "trimmed/{SAMPLES}_rev.fq",
        indexFiles = ["genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    output:
        bams  = "mapped/{SAMPLES}.bam"
    params:
        indexName = "genome/genome"
    message:
        "mapping reads to genome to bam files."
    conda:
        "envs/hisat2.yaml"
    threads: 10
    shell:
        "hisat2 -p {threads} --ignore-quals -x {params.indexName} -1 {input.fwd} -2 {input.rev} | samtools view -Sb -o {output.bams}"

rule create_counts_table:
    input:
        bams = expand("mapped/{sample}.bam", sample=SAMPLES),
        Tct  = "genome/ref_transcriptome.gff"
    output:
        "results/counts.txt"
    #params:
        # some parameters
    conda:
        "envs/subread.yaml"
    shell:
        "featureCounts -F 'gff' -a {input.Tct} -o {output} {input.bams}"

rule DESeq2_analysis:
    input:
        counts    = "results/counts.txt",
        #functions = "results/stringtie_transcriptome_blast.txt"
    output:
        "results/result.csv"
    message:
        "normalizing read counts en creating differential expression table"
    conda:
        "Deseq.yaml"
    shell:
        "Rscript Deseq.R {input.counts}"
