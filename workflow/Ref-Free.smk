# 2024 Benjamin J Perry
# MIT License
# Copyright (c) 2024 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz


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
        expand("results/{library}/06_nonredundant/{samples}.fasta", library = LIBRARY, samples = FIDs), 


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


rule standardise_lengths_nonhost:
    input:
        k2_filtered_reads_gz = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq.gz",
    output:
        filtered_std_reads = "results/{library}/04_k2_filtering/{samples}.nonhost.std.fastq.gz",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "standardise_lengths_nonhost.{samples}.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "standardise_lengths_nonhost.{samples}.txt"),
    conda:
        'cutadapt-4.4'
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 10),
        partition = "compute,hugemem"
    shell:
        """
        cutadapt --cores {threads} --length 65 --minimum-length 65 -o {output.filtered_std_reads} {input.k2_filtered_reads_gz} 
    
        """


rule vsearch_nonredundant_fasta:
    input:
        filtered_std_reads = "results/{library}/04_k2_filtering/{samples}.nonhost.std.fastq.gz",
    output:
        nonredundant_fasta = "results/{library}/04_k2_filtering/{samples}.nonhost.std.nonredundant.fasta",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "vsearch_nonredundant_fasta.{samples}.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "vsearch_nonredundant_fasta.{samples}.txt"),
    conda:
        "vsearch"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 10),
        partition = "compute,hugemem"
    shell:
        "vsearch "
        "--gzip "
        "--threads {threads} "
        "--log {log} "
        "--fastx_uniques {input.filtered_std_reads} "
        "--minuniquesize 1 "
        "--relabel Sequence "
        "--sizeout "
        "--fastaout {output.nonredundant_fasta}"


rule update_nonredundant_headers:
    input:
        nonredundant_fasta = "results/{library}/04_k2_filtering/{samples}.nonhost.std.nonredundant.fasta",
    output:
        nonredundant_fasta_cleaned = "results/{library}/06_nonredundant/{samples}.fasta",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "update_nonredundant_headers.{samples}.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "update_nonredundant_headers.{samples}.txt"),
    threads: 1
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 2 + ((attempt - 1) * 2),
        partition = "compute,hugemem"
    shell:
        """
        
        cat {input.nonredundant_fasta} | sed 's/;/_/g' > {output.nonredundant_fasta_cleaned}
        
        """


# rule report_seqkit_nonredudant:
#     priority: 1000
#     input:
#         expand('results/{library}/01_cutadapt/{samples}.fastq.gz', library = LIBRARY, samples = FIDs),
#     output:
#         os.path.join('results', LIBRARY, '00_QC/seqkit.report.raw.txt')
#     benchmark:
#         os.path.join('benchmarks', LIBRARY, 'report_seqkit_raw.txt')
#     conda:
#         'seqkit'
#     threads: 32
#     resources:
#         mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
#         time = lambda wildcards, attempt: 30 + ((attempt - 1) * 60),
#         partition="compute,hugemem"
#     shell:
#         'seqkit stats -j {threads} -a {input} > {output} '

