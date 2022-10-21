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
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
    #os.system(f"gquery -p no_unpivot -t gbs_keyfile -b library {config['gquery']['libraries']} > resources/{config['gquery']['libraries']}.keyfile.tsv")


#DEPENDENCY - Assumes keyfile is available prior to starting workflow
def getFIDs(keyfile):
    '''
    :param keyfile: gquery generated for library
    :return: list of FIDs
    '''
    import pandas as pd
    keys=pd.read_tsv(keyfile)
    return keys['factid'].tolist()

FIDs = getFIDs('resources/gquery.gbs_keyfile.txt')


rule all:
    input:
        expand('1_cutadapt/{samples}.fastq.gz', samples = FIDs)

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
        '1_cutadapt/{samles}.fastq.gz'
    input: #TODO: determine how to enter the lane level fastq data
        barcodes=rules.generateBarcodes.output.barcodes
    container:
        'docker://...' #TODO
    threads:16
    log:
        'logs/cutadapt.{samples}.fastq.gz'
    params:
    message:
        'Demultiplexing lanes...'
    shell:
        'zcat {input} | '
        'cutadapt '
        '-j 16 '
        '--discard-untrimmed '
        '--no-indels '
        '-g ^file:resources/gquery.barcodes.fasta '
        '-o "1_cutadapt/{name}.fastq.gz" '
        '-' # indicates to use stdin


rule multiQC:
    output:


rule knead
# rule makeFastqLinks:
#     output:



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


