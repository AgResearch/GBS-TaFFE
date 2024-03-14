# 2024 Benjamin J Perry
# MIT License
# Copyright (c) 2024 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

# Based off of Rudi's Homebrew pipeline

configfile: "config/config.yaml"


import os
import pandas as pd


wildcard_constraints:
    samples="\w+"


# Just declaring, in practise this is passed on the CLI by --config LIBRARY=SQ####
LIBRARY = config["LIBRARY"]


input_fastq_pattern = os.path.join('results', config["LIBRARY"], '02_kneaddata', '{samples}.fastq.gz')
print(input_fastq_pattern)
FIDs, = glob_wildcards(input_fastq_pattern)


onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')


rule all:
    input:
        os.path.join("results", LIBRARY, "00_host_stats", (LIBRARY + "_host_multiqc_report.html")),
        os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".host.vcf")),
        os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".merged.host.vcf")),
        os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".homebrew.host.vcf"))


# localrules: get_genome, bcftools_index


rule kraken2_filter:
    input:
        preprocessed_reads = "results/{library}/02_kneaddata/{samples}.fastq.gz",
    output:
        k2OutputHosts = temp("results/{library}/04_k2_filtering/{samples}.filtering.k2"),
        k2_host_reads = temp("results/{library}/04_k2_filtering/{samples}.nonmicrobe.fastq"),
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter.{samples}.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter.{samples}.txt"),
    conda:
        "kraken2"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 420 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
        partition = "hugemem,compute"
    shell:
        "kraken2 "
        "--gzip-compressed "
        "--unclassified-out {output.k2_host_reads} "
        "--db /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1 " 
        "-t {threads} "
        "--output {output.k2OutputHosts} "
        "{input.preprocessed_reads} "
        "2>&1 | tee {log} "


rule kraken2_host_filter_gz:
    input:
        k2_host_reads = "results/{library}/04_k2_filtering/{samples}.nonmicrobe.fastq",
    output:
        k2_host_reads_gz = "results/{library}/04_k2_filtering/{samples}.nonmicrobe.fastq.gz",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:
        "pigz"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 5 + ((attempt - 1) * 5),
        partition = "compute,hugemem"
    shell:
        """

        pigz -p {threads} -c {input.k2_host_reads} > {output.k2_host_reads_gz} &&

        rm {input.k2_host_reads};

        """


# rule get_genome:
#     output:
#         genome_gz = protected('resources/ref/GCF_000298735.2_genomic.fna.gz'),
#     threads: 2
#     resources:
#         partition='compute'
#     resources:
#         time = lambda wildcards, attempt: attempt * 7 * 24 * 60
#     params:
#         genome=config['genome'],
#     shell:
#         "wget -c -O {output.genome_gz} {params.genome} "


# rule build_b2_index:
#     input:
#         genome_gz = 'resources/ref/GCF_000298735.2_genomic.fna.gz',
#     output:
#         bt2_index_semaphore = "resources/ref/BT2INDEX.rev.1.bt2"
#     conda:
#         "bowtie2-2.5.1"
#     threads: 12
#     resources:
#         mem_gb = lambda wildcards, attempt: 16 + ((attempt - 1) * 16), #TODO benchmark and tweak
#         time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
#         partition = "compute"
#     shell:
#         """
        
#         bowtie2 -f --threads {threads} --seed 1953 {input.genome_gz} resources/ref/BT2INDEX 
#         bowtie2-inspect -s resources/ref/BT2INDEX

#         """


rule bowtie2_alignment:
    input:
        k2_host_reads_gz = "results/{library}/04_k2_filtering/{samples}.nonmicrobe.fastq.gz",
        bt2_index_semaphore = "resources/ref/BT2INDEX.rev.1.bt2"
    output:
        host_sam = temp("results/{library}/06_host_alignment/{samples}.sam"),
    log:
        os.path.join("results", "{library}", "logs", "bowtie2", "bowtie2_alignment.{samples}.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "bowtie2_alignment.{samples}.txt"),
    conda:
        "bowtie2-2.5.1"
    threads: 6
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 84),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 20),
        partition = "compute,hugemem"
    shell:
        "bowtie2 "
        "--very-fast-local "
        "-p {threads} "
        "--rg-id {wildcards.samples} "
        "--rg SM:{wildcards.samples} "
        "--rg LB:{wildcards.library} "
        "--rg PU:{wildcards.library} " #TODO add flowcell here
        "--rg PL:ILLUMINA " 
        "-x resources/ref/BT2INDEX "
        "-U {input.k2_host_reads_gz} "
        "-S {output.host_sam} "
        "2> {log} "


rule prepare_bams:
    input:
        host_sam = "results/{library}/06_host_alignment/{samples}.sam",
    output:
        host_bam = "results/{library}/06_host_alignment/{samples}.sorted.bam",
    log:
        os.path.join("results", "{library}", "logs", "samtools", "prepare_bams.{samples}.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "prepare_bams.{samples}.txt"),
    conda:
        "samtools-1.17"
    threads: 6
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 60),
        partition = "compute,hugemem"
    shell:
        """
        samtools view --threads 6 -bS  {input.host_sam} | samtools sort > {output.host_bam} 2> {log}


        """


# rule bamtools_merge_bams:
#     input:
#         host_bams = expand("results/{library}/06_host_alignment/{samples}.sorted.bam", library = LIBRARY, samples = FIDs),
#     output:
#         bam_list = os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".bamlist.txt")),
#         merged_bams = os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".merged.host.bamtools.bam"))
#     log:
#         os.path.join("results", LIBRARY, "logs", "bamtools", "bcftools_merge_bams.log"),
#     benchmark:
#         os.path.join("results", LIBRARY, "benchmarks", "bcftools_merge_bams.txt"),
#     conda:
#         "bamtools-2.5.2"
#     threads: 12
#     resources:
#         mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 32),
#         time = lambda wildcards, attempt: 720 + ((attempt - 1) * 720),
#         partition = "compute"
#     shell:
#         """
#         for i in {input.host_bams}; do echo $i >> {output.bam_list}; done &&
        
#         bamtools merge -list {output.bam_list} -out {output.merged_bams}

#         """


rule samtools_merge_bams:
    input:
        host_bams = expand("results/{library}/06_host_alignment/{samples}.sorted.bam", library = LIBRARY, samples = FIDs),
        bcf_index = 'resources/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna' #TODO automate the file name expansion to add .fai
    output:
        merged_bams = os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".merged.host.samtools.bam"))
    log:
        os.path.join("results", LIBRARY, "logs", "bamtools", "samtools_merge_bams.log"),
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "samtools_merge_bams.txt"),
    conda:
        "samtools-1.17"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 720 + ((attempt - 1) * 720),
        partition = "compute"
    shell:
        "samtools merge "
        "-o {output.merged_bams} "
        "-u "
        "-O BAM "
        "--reference {input.bcf_index} "
        "--threads {threads} "
        "--write-index "
        "{input.host_bams} "


# rule bcftools_index:
#     input:
#         genome_gz = 'resources/ref/GCF_000298735.2_genomic.fna.gz',
#     output:
#         bcf_index = 'resources/ref/GCF_000298735.2_genomic.fna' #TODO automate the file name expansion to add .fai
#     conda:
#         "samtools-1.17"
#     threads: 1
#     shell:
#         """
#         wget -k {input.genome_gz};
#         samtools faidx resources/ref/GCF_000298735.2_genomic.fna
        
#         """


rule bcftools_VCF_merged: #Alternate Method for VCF
    input:
        merged_bams = os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".merged.host.samtools.bam")),
        indexed_genome = 'resources/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna',
    output:
        host_vcf = os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".merged.host.vcf")),
    log:
        os.path.join("results", LIBRARY, "logs", "bcftools", "bcftools_VCF.log"),
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "bcftools_VCF.txt"),
    conda:
        "bcftools-1.19"
    threads: 24
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 2880 + ((attempt - 1) * 720),
        partition = "compute"
    shell:
        "bcftools mpileup --threads {threads} -I -Ou -f {input.indexed_genome} -a INFO/DPR,INFO/AD,FORMAT/DP,FORMAT/AD {input.merged_bams} "
        "| bcftools call -cv - "
        "| bcftools view -M2 - "
        "> {output.host_vcf} "


rule bcftools_VCF_homebrew: #Alternate Method for VCF
    input:
        host_bams = expand("results/{library}/06_host_alignment/{samples}.sorted.bam", library = LIBRARY, samples = FIDs),
        indexed_genome = 'resources/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna',
    output:
        bam_list = os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".bamlist.txt")),
        host_vcf = os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".homebrew.host.vcf")),
    log:
        os.path.join("results", LIBRARY, "logs", "bcftools", "bcftools_VCF_list.log"),
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "bcftools_VCF_list.txt"),
    conda:
        "bcftools-1.19"
    threads: 24
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 12),
        time = lambda wildcards, attempt: 2880 + ((attempt - 1) * 720),
        partition = "compute"
    shell:
        """
        for i in {input.host_bams}; do echo $i >> {output.bam_list}; done && 
        
        bcftools mpileup -I -Ou -f {input.indexed_genome} -b {output.bam_list} -a AD | bcftools call -cv - | bcftools view -M2 - > {output.host_vcf};

        """


rule bcftools_VCF_individual:
    input:
        host_bam = "results/{library}/06_host_alignment/{samples}.sorted.bam",
        indexed_genome = 'resources/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna',
    output:
        host_vcf = "results/{library}/06_host_alignment/{samples}.sorted.bam.vcf.gz",
        csi = "results/{library}/06_host_alignment/{samples}.sorted.bam.vcf.gz.csi",
    log:
        "results/{library}/logs/bcftools/bcftools_VCF_individual.{samples}.log",
    benchmark:
        "results/{library}/benchmarks/bcftools_VCF_individual.{samples}.txt",
    conda:
        "bcftools-1.19"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 6 + ((attempt - 1) * 24),
        partition = "compute"
    shell:
        "bcftools mpileup --threads {threads} -I -Ou -f {input.indexed_genome} -a INFO/DPR,INFO/AD,FORMAT/DP,FORMAT/AD {input.host_bam} "
        "| bcftools call -cv - "
        "| bcftools view -M2 -O z8 -o {output.host_vcf}; "
        "bcftools index --threads {threads} {output.host_vcf} -o  {output.csi} "


rule merge_bcftools_VCF_individual:
    input:
        vcf = expand("results/{library}/06_host_alignment/{samples}.sorted.bam.vcf.gz", samples = FIDs, library = LIBRARY),
        csi = expand("results/{library}/06_host_alignment/{samples}.sorted.bam.vcf.gz.csi", samples = FIDs, library = LIBRARY),
    output:
        host_vcf = os.path.join("results", LIBRARY, "06_host_alignment", (LIBRARY + ".host.vcf")),
    log:
        os.path.join("results", LIBRARY, "logs", "bcftools", "merge_bcftools_VCF_individual.log"),
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "merge_bcftools_VCF_individual.txt"),
    conda:
        "bcftools-1.19"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 12),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute"
    shell:
        """
        bcftools merge --threads {threads} {input.vcf} | bcftools view -m2 -M2 -v snps -O v -o {output.host_vcf};
        """


rule samtools_stats:
    input:
        host_bam = "results/{library}/06_host_alignment/{samples}.sorted.bam",
        reference = 'resources/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna' #TODO automate the file name expansion
    output:
        stats = "results/{library}/00_host_stats/{samples}.bam.samtools_stats.txt",
    log:
        "results/{library}/logs/bcftools/samtools_stats.{samples}.log",
    benchmark:
        "results/{library}/benchmarks/samtools_stats.{samples}.txt",
    conda:
        "samtools-1.17"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 4 + ((attempt - 1) * 60),
        partition = "compute",
    shell:
        "samtools stats --threads {threads} -r {input.reference} {input.host_bam} > {output.stats} "


rule samtools_index_bam:
    input:
        host_bam = "results/{library}/06_host_alignment/{samples}.sorted.bam",
    output:
        index = "results/{library}/06_host_alignment/{samples}.sorted.bam.bai",
    log:
        "results/{library}/logs/bcftools/{samples}.samtools_index_bam.log",
    benchmark:
        "results/{library}/benchmarks/{samples}.samtools_index_bam.txt",
    conda:
        "samtools-1.17"
    threads: 1
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 2 + ((attempt - 1) * 10),
        partition = "compute",
    shell:
        """
        samtools index {input.host_bam};
        """


rule mosdepth_stats:
    input:
        host_bam = "results/{library}/06_host_alignment/{samples}.sorted.bam",
        index = "results/{library}/06_host_alignment/{samples}.sorted.bam.bai",
    output:
        stats = "results/{library}/00_host_stats/{samples}.host.mosdepth.summary.txt",
    log:
        "results/{library}/logs/bcftools/mosdepth_stats.{samples}.log",
    benchmark:
        "results/{library}/benchmarks/mosdepth_stats.{samples}.txt",
    conda:
        "mosdepth-0.3.6"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 2 + ((attempt - 1) * 2),
        partition = "compute",
    shell:
        "mosdepth "
        "--use-median "
        "--fast-mode " # dont look at internal cigar operations or correct mate overlaps (recommended for most use-cases).
        "--no-per-base " # dont output per-base depth.
        "--threads {threads} "
        "results/{wildcards.library}/00_host_stats/{wildcards.samples}.host " # output prefix
        "{input.host_bam} "


rule bcftools_stats:
    input:
        host_vcf = "results/{library}/06_host_alignment/{samples}.sorted.bam.vcf.gz",
        reference = 'resources/ref/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna' #TODO automate the file name expansion
    output:
        stats = "results/{library}/00_host_stats/{samples}.host.bcftools-stats.txt",
    log:
        "results/{library}/logs/bcftools/bcftools_stats.{samples}.log",
    benchmark:
        "results/{library}/benchmarks/bcftools_stats.{samples}.txt",
    threads: 4
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 2),
        time = lambda wildcards, attempt: 2 + ((attempt - 1) * 2),
        partition = "compute",
    shell:
        "bcftools stats "
        "--fasta-ref {input.reference} "
        "--samples - "
        "--threads {threads} "
        "{input.host_vcf} > "
        "{output.stats} "


rule host_multiqc:
    input:
        logs_bowtie2 = expand("results/{library}/logs/bowtie2/bowtie2_alignment.{samples}.log", library = LIBRARY, samples = FIDs),
        stats_bcftools = expand("results/{library}/00_host_stats/{samples}.host.bcftools-stats.txt", library = LIBRARY, samples = FIDs),
        stats_mosdepth = expand("results/{library}/00_host_stats/{samples}.host.mosdepth.summary.txt", library = LIBRARY, samples = FIDs),
        stats_samtools = expand("results/{library}/00_host_stats/{samples}.bam.samtools_stats.txt", library = LIBRARY, samples = FIDs),
    output:
        multiqc_report = os.path.join("results", LIBRARY, "00_host_stats", (LIBRARY + "_host_multiqc_report.html")),
    log:
        os.path.join("results", LIBRARY, "logs", "multiqc", "host_multiqc.log"),
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "host_multiqc.txt"),
    conda:
        "multiqc"
    threads: 6
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
        partition = "compute",
    shell:
        "multiqc "
        "--interactive "
        "--title {LIBRARY}_host "
        "--force "
        "--data-format tsv "
        "--fullnames "
        "--outdir results/{LIBRARY}/00_host_stats "
        "{input.logs_bowtie2} {input.stats_bcftools} {input.stats_mosdepth} {input.stats_samtools} "

