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
if not os.path.isdir(library_keyfile):
    print("Generating Keyfile for " + library_keyfile)
    shell("mkdir -p resources")
    shell(f"gquery -p no_unpivot -t gbs_keyfile -b library {config['gquery']['libraries']} > {library_keyfile}")


#DEPENDENCY - Assumes keyfile is available prior to starting workflow
def getFIDs(keyfile):
    '''
    :param keyfile: gquery generated for library
    :return: list of FIDs
    '''
    import pandas as pd
    keys=pd.read_tsv(keyfile)
    return keys['factid'].tolist()

FIDs = getFIDs(library_keyfile)
print(FIDs)

rule all:
    input:
        expand('01_cutadapt/{samples}.fastq.gz', samples = FIDs),
        '0_qc/mergedReadsMultiQCReport.html'

### Currently Assuming Keyfile is Already: resources/gquery.keyfile.txt
# rule generateKeyfile: #TODO replace with rule generateBarcodes when gquery has feature
#     output:
#         keyfile = 'resources/gquery.gbs_keyfile.txt'
#     threads:2
#     log:
#         'logs/1_gqueryGenerateKeyfile.log'
#     params:
#         libraries = config['gquery']['libraries'],
#     message:
#         'Dumping keyfile for: {params.libraries}...\n'
#     shell:
#         'gquery '
#         '-p no_unpivot '
#         '-t gbs_keyfile '
#         '-b library {params.libraries} > '
#         '{output.keyfile} 2> '
#         '{log}'


rule generateBarcodes: #TODO replace with rule generateBarcodes when gquery has feature
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
        expand('01_cutadapt/{samples}.fastq.gz', samples = FIDs)
    input:
        barcodes = rules.generateBarcodes.output.barcodes,
        lane01 = config['novaseq']['lane01'],
        lane02 = config['novaseq']['lane02'],
    container:
        'docker://...' #TODO
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
        r'-o "1_cutadapt/{{name}}.fastq.gz" '
        '-' # indicates stdin


rule fastqcMerged:
    output:
        html = '0_qc/fastqc/{samples}_fastqc.html',
        zip = "0_qc/fastqc/{samples}_fastqc.zip"
    input:
        fastq = '01_cutadapt/{samples}.fastq.gz'
    conda:
        'fastqc'
    threads: 4
    message:
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o 0_qc/fastqc/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule multiQCMerged:
    output:
        multiQC='0_qc/mergedReadsMultiQCReport.html'
    input:
        fastqc= expand('0_qc/fastqc/{samples}_fastqc.zip', samples = SAMPLES)
    conda:
        'multiqc'
    shell:
        'multiqc '
        '-n 0_qc/mergedReadsMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'


# rule kneaddata:


# rule humann3:


# rule human3GTDB:


# rule metaphlan4:


# rule kraken2:


# rule kraken2GTDB:



# ### Exapmle Rule Template ###
# rule cutadapt:
#     input:
#         'fastq/{samples}.fastq.gz'
#     output:
#         '02_cutadapt/{samples}.chop.primer.fastq.gz'
#     threads:4
#     log:
#         'logs/{samples}/cutadapt.primers.log'
#     container: 
#         'docker://quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'
#     params:
#         fwdPrimer=config['cutadapt']['fwd'],
#         revPrimer=config['cutadapt']['rev'],
#     resources:
#         tempdir=config['TMPDIR']
#     message:
#         'removing primers: {wildcards.samples}\n'
#         'TMPDIR: {resources.tempdir}'
#     shell:
#         'cutadapt '
#         '--discard-untrimmed '
#         '--action=retain '
#         '-j {threads} '
#         '--error-rate 0.2 '
#         '-g {params.fwdPrimer} '
#         '-a {params.revPrimer} '
#         '-o {output} '
#         '{input} '
#         '2>&1 | tee {log}'


