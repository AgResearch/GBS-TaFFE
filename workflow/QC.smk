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


FIDs, = glob_wildcards("fastq/{sample}_L001_R1_001.fastq.gz") #TODO


rule all:
    input:
        'results/00_QC/MultiQCReport.html',


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

