# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"

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


library_keyfile=f"resources/{config['gquery']['libraries']}.keyfile.tsv"
# if no keyfile present in resources/, spew keyfile for GBS library
if not os.path.isfile(library_keyfile):
    print("Generating Keyfile: " + library_keyfile)
    shell("mkdir -p resources")
    shell(f"gquery -p no_unpivot -t gbs_keyfile -b library {config['gquery']['libraries']} > {library_keyfile}")


# Extracthe unique factids for samples to use as 'wildcards'
def getFIDs(keyfile):
    '''
    :param keyfile: gquery generated for library
    :return: list of FIDs
    '''
    import pandas as pd
    keys=pd.read_csv(keyfile, sep='\t', header=0)
    return keys['factid'].tolist()

print(f"Extracting factid for {config['gquery']['libraries']}...")
FIDs = getFIDs(library_keyfile)

print("Found: ")
for entry in FIDs:
    print(entry)


rule all:
    input:
        expand('02_kneaddata/{samples}_kneaddata.fastq.gz', samples = FIDs),
        '00_qc/ReadsMultiQCReport.html',
        '00_qc/KDRReadsMultiQCReport.html'


rule generateBarcodes:
    output:
        barcodes = 'resources/gquery.barcodes.fasta'
    threads:2
    log:
        'logs/1_gqueryGenerateBarcodes.log'
    params:
        libraries = config['gquery']['libraries'],
    message:
        'Spewing barcodes for: {params.libraries}...\n'
    shell:
        'gquery '
        '-t gbs_keyfile '
        '-b library -p "columns=factid,barcode;fasta;noheading;no_unpivot" '
        '{params.libraries} > '
        '{output.barcodes} 2> '
        '{log}'


rule cutadapt:
    output:
        expand('01_cutadapt/{samples}.fastq.gz', samples = FIDs),
    input:
        barcodes = rules.generateBarcodes.output.barcodes,
        lane01 = config['novaseq']['lane01'],
        lane02 = config['novaseq']['lane02'],
    container:
        'docker://quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'
    threads:16
    log:
        'logs/cutadapt.log'
    params:
    message:
        'Demultiplexing lanes...'
    shell:
        'zcat {input.lane01} {input.lane02} | '
        'cutadapt '
        '-j 16 '
        '--discard-untrimmed '
        '--no-indels '
        '-g ^file:{input.barcodes} '
        r'-o "01_cutadapt/{{name}}.fastq.gz" '
        '-' # indicates stdin


rule fastqc:
    output:
        html = '00_qc/fastqc/{samples}_fastqc.html',
        zip = "00_qc/fastqc/{samples}_fastqc.zip"
    input:
        fastq = '01_cutadapt/{samples}.fastq.gz'
    container:
        'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 1
    message:
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o 00_qc/fastqc/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule multiQC:
    output:
        multiQC='00_qc/ReadsMultiQCReport.html'
    input:
        fastqc= expand('00_qc/fastqc/{samples}_fastqc.zip', samples = FIDs)
    container:
        'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '-n 00_qc/ReadsMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'

#TODO Rule to Build Rambv2 index

rule kneaddata:
    input:
        reads = '01_cutadapt/{samples}.fastq.gz',
    output:
        trimReads = temp('02_kneaddata/{samples}_kneaddata.trimmed.fastq'),
        trfReads = temp('02_kneaddata/{samples}_kneaddata.repeats.removed.fastq'),
        ovineReads = temp('02_kneaddata/{samples}_kneaddata_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq'),
        clnReads= '02_kneaddata/{samples}_kneaddata.fastq',
        readStats = '02_kneaddata/{samples}.read.stats.txt'
    log:
        'logs/{samples}.kneaddata.log'
    conda:
        'biobakery'
    threads: 4
    message:
        'kneaddata: {wildcards.samples}\n'
    shell:
        'kneaddata '
        '--input {input.reads} '
        '-t {threads} '
        '--log-level INFO '
        '--log {log} '
        '--trimmomatic /home/perrybe/conda-envs/biobakery/share/trimmomatic '
        '--sequencer-source TruSeq3 '
        '-db ref/Rambv2/GCF_016772045.1-ARS-UI-Ramb-v2.0 '
        '-o 02_kneaddata && '
        'seqkit stats -j 12 -a 02_kneaddata/{wildcards.samples}*.fastq > {output.readStats}'

#TODO Compress output reads

rule fastqcKDRs:
    output:
        html = '00_qc/fastqcKDR/{samples}_fastqc.html',
        zip = "00_qc/fastqcKDR/{samples}_fastqc.zip"
    input:
        fastq = rules.kneaddata.output.clnReads
    container:
        'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 1
    message:
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o 00_qc/fastqcKDR/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule multiQCKDRs:
    output:
        multiQC='00_qc/KDRReadsMultiQCReport.html'
    input:
        fastqc= expand('00_qc/fastqcKDR/{samples}_fastqc.html', samples = FIDs)
    container:
        'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '-n 00_qc/KDRReadsMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'

# rule kraken2:
#     output:

#     input:

#     log:

#     threads:

#     message:

#     shell:

# rule kraken2GTDB:
#     output:

#     input:

#     log:

#     threads:

#     message:

#     shell:


# rule humann3:
#     output:

#     input:

#     log:

#     threads:

#     message:

#     shell:


# rule human3GTDB:
#     output:

#     input:

#     log:

#     threads:

#     message:

#     shell:


# rule metaphlan4:
#     output:

#     input:

#     log:

#     threads:

#     message:

#     shell:

