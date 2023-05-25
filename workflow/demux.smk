# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

import os

LIBRARIES, = glob_wildcards("resources/{library}.barcodes.fasta")

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
        expand("results/{library}", library=LIBRARIES),

rule cutadapt: # demultiplexing GBS reads
    input:
        barcodes = "resources/{library}.barcodes.fasta",
        run = "fastq/{library}.fastq.gz",
    output:
        demuxed = directory("results/{library}"),
    container:
        'docker://quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'
    benchmark:
        'benchmarks/cutadapt.{library}.txt'
    threads: 32
    resources:
        mem_gb=8,
        time="01:00:00",
	partition="large,milan"
    message:
        'Demultiplexing lanes...'
    shell:
        'mkdir -p {output.demuxed} && '
        'zcat {input.run} | '
        'cutadapt '
        '-j {threads} '
        '--discard-untrimmed '
        '--no-indels '
        '-g ^file:{input.barcodes} '
        r'-o "{output.demuxed}/{{name}}.fastq.gz" '
        '-'

