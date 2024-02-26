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
        expand("results/{library}/{library}.kraken2.GTDB214.domain.counts.tsv",  library = LIBRARY),


rule kraken2_host_filter:
    input:
        preprocessed_reads = "results/{library}/02_kneaddata/{samples}.fastq.gz",
    output:
        k2OutputHosts = temp("results/{library}/04_k2_filtering/{samples}.hosts.k2"),
        k2_filtered_reads = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq",
        k2_host_reads = "results/{library}/04_k2_filtering/{samples}.host.fastq",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter.{samples}.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter.{samples}.txt"),
    conda:
        "kraken2"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 32), #TODO benchmark and tweak
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
        partition = "compute,hugemem"
    shell:
        "kraken2 "
        "--gzip-compressed "
        "--unclassified-out {output.k2_filtered_reads} "
        "--classified-out {output.k2_host_reads} "
        "--db /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-hosts " 
        "-t {threads} "
        "--output {output.k2OutputHosts} "
        "{input.preprocessed_reads} "
        "2>&1 | tee {log} "


rule kraken2_host_filter_gz:
    input:
        k2_filtered_reads = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq",
        k2_host_reads = "results/{library}/04_k2_filtering/{samples}.host.fastq",
    output:
        k2_filtered_reads_gz = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq.gz",
        k2_host_reads_gz = "results/{library}/04_k2_filtering/{samples}.host.fastq.gz",
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

        pigz -p {threads} -c {input.k2_filtered_reads} > {output.k2_filtered_reads_gz} && 
        pigz -p {threads} -c {input.k2_host_reads} > {output.k2_host_reads_gz} &&

        rm {input.k2_filtered_reads} {input.k2_host_reads};

        """


rule bowtie2_alignment:
    input:

    output:

    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:

    threads:

    shell:
        """
        
        """


rule prepare_bams:
    input:

    output:

    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:

    threads:

    shell:
        """
        
        """


rule bcftools_VCF:
    input:

    output:

    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:

    threads:

    shell:
        """
        
        """


rule normalise_VCF:
    input:

    output:

    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:

    threads:

    shell:
        """
        
        """


rule merge_VCF:
    input:

    output:

    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:

    threads:

    shell:
        """
        
        """


