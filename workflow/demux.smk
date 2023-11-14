# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

import os


LIBRARIES, = glob_wildcards("resources/{library}.cutadapt.barcodes.fasta")


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
        expand("results/01_cutadapt/{library}", library=LIBRARIES),


rule cutadapt: # demultiplexing GBS reads
    input:
        barcodes = "resources/{library}.cutadapt.barcodes.fasta",
    output:
        demuxed = directory("results/01_cutadapt/{library}"),
    conda:
        'cutadapt-4.4'
    benchmark:
        'benchmarks/cutadapt.{library}.txt'
    threads: 32
    resources:
        mem_gb=24,
        time=240,
	    partition="compute,hugemem"
    shell:
        'mkdir -p {output.demuxed} && '
        'zcat fastq/{wildcards.library}*.gz | '
        'cutadapt '
        '-j {threads} '
        '--discard-untrimmed '
        '--no-indels '
        '-e 0 '
        '-g ^file:{input.barcodes} '
        r'-o "{output.demuxed}/{{name}}.fastq.gz" '
        '-'

