configfile: "config.yaml"

# fetch URL to transcriptome multi fasta from configfile
genome_URL = config["refs"]["genome_URL"]
transcriptome_URL = config["refs"]["transcriptome_URL"]
# create lists containing samplenames, conditions, path/filenames of the fw-reads(R1)
# and the path/filenames of the rev reads(R2) from the file: data/sampls.txt
import pandas as pd
samples = list(pd.read_table("data/samples.txt")["sample"])
conditions = list(pd.read_table("data/samples.txt")["condition"])
R1 = list(pd.read_table("data/samples.txt")["fq1"])
R2 = list(pd.read_table("data/samples.txt")["fq2"])

rule all:
    output:
        "results/result.csv"
        fwd = "fastqc/{samples}_forward_fastqc.zip",
        rev = "fastqc/{samples}_reverse_fastqc.zip"
    massage:
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
    shell: "wget -O {output} {transcriptome_URL}"

rule get_ref_transcriptome_index:
    input:
        "genome/ref_transcriptome.fasta"
    output:
        ["genome/ref_transcriptome.fasta." + i for i in ("nhr", "nin", "nog", "nsd", "nsi", "nsd")]
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

rule index:
    input:
        "genome/genome.fasta"
    output:
        ["genome/genome." + str(i) + ".bt2" for i in range(1,5)],
        "genome/genome.rev.1.bt2",
        "genome/genome.rev.2.bt2"
    message:"indexing genome"
    params:
        "genome/genome"
    threads: 10
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule trimmomatic:
    input:
        fq1 = R1,
        fq2 = R2,
        adapters = config["adapters"]
    output:
        fw_reads = "trimmed/{sample}_fw.fq.gz",
        rev_reads = "trimmed/{sample}_rev.fq.gz",
        forwardUnpaired = temp("trimmed/{sample}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp("trimmed/{sample}_reverse_unpaired.fastq.gz")
    message: "trimming {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic/{sample}_{unit}.log"
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
        "{input.reads} "
        "{output.forward_reads} "
        "{output.forwardUnpaired} "
        "{output.reverse_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} 2>{log}"

rule fastqc:
    input:
        fwd = "trimmed/{sample}_fw.fq.gz",
        rev = "trimmed/{sample}_rev.fq.gz",
    output:
        fwd="fastqc/{samples}_forward_fastqc.zip",
        rev="fastqc/{samples}_reverse_fastqc.zip"
    log:
        "results/logs/fastqc/{samples}.fastqc.log"
    params:
        "results/fastqc/"
    message:
        "Quality check of trimmed samples with FASTQC" 		#removed, it was not working
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev} 2>{log}"

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
