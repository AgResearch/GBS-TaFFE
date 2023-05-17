# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

#configfile: "config/config.yaml"


import os


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


FIDs, = glob_wildcards("results/03_kraken2GTDB/{sample}.kraken2") #TODO


rule all:
    input:
        'results/00_QC/MultiQCReport.html',


rule fastqc:
    input:
        fastq = 'results/01_cutadapt/{samples}.fastq.gz'
    output:
        html = 'results/00_QC/fastqc/{samples}_fastqc.html',
        zip = 'results/00_QC/fastqc/{samples}_fastqc.zip'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    benchmark:
        'benchmarks/fastqc.{samples}.txt'
    threads: 2
    message:
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o results/00_QC/fastqc/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule fastqcKDRs:
    input:
        fastq = 'results/02_kneaddata/{samples}.fastq'
    output:
        'results/00_QC/fastqcKDR/{samples}_fastqc.zip'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    benchmark:
        'benchmarks/fastqcKDRs.{samples}.txt'
    threads: 2
    message:
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o results/00_QC/fastqcKDR/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule multiQC:
    input:
        fastqc = expand('results/00_QC/fastqc/{sample}_fastqc.zip', sample = FIDs),
        KDRs = expand('results/00_QC/fastqcKDR/{sample}_fastqc.zip', sample = FIDs),
        bbduk = expand('logs/bbduk/{sample}.bbduk.log', sample = FIDs),
        prinseq = expand('logs/prinseq/{sample}.prinseq.log', sample = FIDs),
        kneaddata = expand('logs/kneaddata/{sample}.kneaddata.log', sample = FIDs),
        kraken2 = expand('results/03_kraken2GTDB/{sample}.kraken2', sample = FIDs)
    output:
        multiQC ='results/00_QC/MultiQCReport.html'
    conda:
        'multiqc'
        # 'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '-n results/00_QC/MultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc} {input.KDRs} {input.bbduk} {input.prinseq} {input.kneaddata} {input.kraken2} '

