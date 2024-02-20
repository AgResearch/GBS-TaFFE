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


# Global minimum read count for processing
min_reads = 25000

LIBRARY = config["LIBRARY"]

# seqkit_report = os.path.join("results", LIBRARY, "00_QC", "seqkit.report.raw.txt")

# def get_passing_FIDs(seqkitOut = seqkit_report, minReads=min_reads, lib=LIBRARY):
#     import pandas as pd
#     qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
#     qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
#     qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
#     return qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()

# FIDs = get_passing_FIDs()


input_fastq_pattern = os.path.join('results', config["LIBRARY"], '01_cutadapt', '{samples}.fastq.gz')
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
        k2OutputHosts = temp("results/{library}/04_k2_filtering/{samples}.Hosts.k2"),
        k2_filtered_read = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq",
        k2_host_read = "results/{library}/04_k2_filtering/{samples}.host.fastq",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter.{samples}.txt"),
    conda:
        "kraken2"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
        partition = "compute,hugemem"
    shell:
        "kraken2 "
        "--gzip-compressed "
        "--unclassified-out {output.k2_filtered_read} "
        "--classified-out {output.k2Classified_read} "
        "--db /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-hosts " 
        "-t {threads} "
        "--output {output.k2OutputHosts} "
        "{input.k2Classified_read} "
        "2>&1 | tee {log} "


rule kraken2_host_filter_gz:
    input:
        k2_filtered_read = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq",
        k2_host_read = "results/{library}/04_k2_filtering/{samples}.host.fastq",
    output:
        k2_filtered_reads_gz = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq.gz",
        k2_host_read_gz = "results/{library}/04_k2_filtering/{samples}.host.fastq.gz",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:
        "pigz"
    threads: 16
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 16),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 20),
        partition = "compute,hugemem"
    shell:
        """

        pigz -p {threads} -c {input.k2_filtered_read} > {output.k2_filtered_reads_gz} && 
        pigz -p {threads} -c {input.k2_host_read} > {output.k2_host_read_gz} &&

        rm {input.k2_filtered_read} {input.k2_host_read};

        """


rule standardise_length:
    input:
        k2_filtered_reads_gz = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq.gz",
    output:
        filtered_std_reads = "results/{library}/04_k2_filtering/{samples}.nonhost.std.fastq.gz",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:
        'cutadapt-4.4'
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 16),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 20),
        partition = "compute,hugemem"
    shell:
        """

    
        """


rule vsearch_dereplicate:
    input:

    output:
        nonredundant_fasta = "results/{library}/05_nonredundant/{samples}.nonhost.fasta",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:
        "vsearch"
    threads: 1
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 16),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 20),
        partition = "compute,hugemem"
    shell:
        "vsearch "
        "--gzip "
        "--fastaout "
        "--threads {threads} "
        "--log {log} "
        "--fastx_uniques {input.k2_filtered_reads_gz} "
        "--minuniquesize 1 "
        "--relabel Sequence "
        "--sizeout "
        "--fastaout {output.nonredundant_fasta}"


rule report_seqkit:
    priority: 1000
    input:
        expand('results/{library}/01_cutadapt/{samples}.fastq.gz', library = LIBRARY, samples = FIDs),
    output:
        os.path.join('results', LIBRARY, '00_QC/seqkit.report.raw.txt')
    benchmark:
        os.path.join('benchmarks', LIBRARY, 'report_seqkit_raw.txt')
    conda:
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 60),
        partition="compute,hugemem"
    shell:
        'seqkit stats -j {threads} -a {input} > {output} '

