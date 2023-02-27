# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"

import os

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
    print("Found: ")
    for entry in FIDs:
        print(entry)


	
rule all:
    input:
        'results/00_qc/ReadsMultiQCReport.html',
        'results/00_qc/KDRReadsMultiQCReport.html',
        expand('results/02_kneaddata/{samples}_kneaddata.trimmed.fastq', samples=FIDs),
        'results/00_qc/seqkit.report.raw.txt',
        'results/00_qc/seqkit.report.filtered.txt',




localrules: generateBarcodes



rule generateBarcodes:
    output:
        barcodes = 'resources/gquery.barcodes.fasta'
    threads: 2
    log:
        'logs/1_gqueryGenerateBarcodes.log'
    params:
        libraries = config['gquery']['libraries'],
    message:
        'Generating barcode file for: {params.libraries}...\n'
    shell:
        'gquery '
        '-t gbs_keyfile '
        '-b library -p "columns=factid,barcode;fasta;noheading;no_unpivot" '
        '{params.libraries} > '
        '{output.barcodes} 2> '
        '{log}'



rule cutadapt: # demultiplexing GBS reads
    input:
        barcodes = rules.generateBarcodes.output.barcodes,
        lane01 = config['novaseq']['lane01'],
        lane02 = config['novaseq']['lane02'],
    output:
        expand('results/01_cutadapt/{samples}.fastq.gz', samples = FIDs),
    conda:
        'cutadapt'
        # 'docker://quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'
    threads: 32
    resources:
        mem_gb=64,
        partition='inv-iranui-fast'
    message:
        'Demultiplexing lanes...'
    shell:
        'zcat {input.lane01} {input.lane02} | '
        'cutadapt '
        '-j {threads} '
        '--discard-untrimmed '
        '--no-indels '
        '-g ^file:{input.barcodes} '
        r'-o "results/01_cutadapt/{{name}}.fastq.gz" '
        '-' # indicates stdin



rule fastqc:
    input:
        fastq = 'results/01_cutadapt/{samples}.fastq.gz'
    output:
        html = 'results/00_qc/fastqc/{samples}_fastqc.html',
        zip = 'results/00_qc/fastqc/{samples}_fastqc.zip'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
    message:
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o results/00_qc/fastqc/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'



rule multiQC:
    input:
        fastqc= expand('results/00_qc/fastqc/{samples}_fastqc.zip', samples = FIDs)
    output:
        multiQC='results/00_qc/ReadsMultiQCReport.html'
    conda:
        'multiqc'
        # 'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '-n results/00_qc/ReadsMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'



rule kneaddata:
    input:
        reads = 'results/01_cutadapt/{samples}.fastq.gz',
    output:
        #trimReads = temp('results/02_kneaddata/{samples}_kneaddata.trimmed.fastq'),
        #trfReads = temp('results/02_kneaddata/{samples}_kneaddata.repeats.removed.fastq'),
        #ovineReads = temp('results/02_kneaddata/{samples}_kneaddata_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq'),
        #silvaReads = temp('results/02_kneaddata/{samples}_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam.fastq'),
        KDRs = 'results/02_kneaddata/{samples}_kneaddata.trimmed.fastq',
    conda:
        'biobakery'
    log:
        'logs/kneaddata/{samples}.kneaddata.log'
    threads: 8
    resources:
        mem_gb=8,
        time='02:00:00'
    message:
        'kneaddata: {wildcards.samples}\n'
    shell:
        'kneaddata '
        '--trimmomatic-options "ILLUMINACLIP:/home/perrybe/conda-envs/biobakery/share/trimmomatic-0.39-2/adapters/illuminaAdapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:40" '
        '--input {input.reads} '
        '-t {threads} '
        '--log-level INFO '
        '--log {log} '
        '--trimmomatic /home/perrybe/conda-envs/biobakery/share/trimmomatic '
        '--sequencer-source TruSeq3 '
        '--bypass-trf '
        #'-db /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/Rambv2/GCF_016772045.1-ARS-UI-Ramb-v2.0 '
        #'-db /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/SILVA_128_LSUParc_SSUParc_ribosomal_RNA '
        '-o results/02_kneaddata '



rule fastqcKDRs:
    input:
        fastq = 'results/02_kneaddata/{samples}_kneaddata.trimmed.fastq'
    output:
        'results/00_qc/fastqcKDR/{samples}_kneaddata.trimmed_fastqc.zip'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
    message:
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o results/00_qc/fastqcKDR/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'



rule multiQCKDRs:
    input:
        fastqc= expand('results/00_qc/fastqcKDR/{samples}_kneaddata.trimmed_fastqc.zip', samples = FIDs)
    output:
        'results/00_qc/KDRReadsMultiQCReport.html'
    conda:
        'multiqc'
        # 'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '-n results/00_qc/KDRReadsMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'



rule seqkitQCRaw:
    input:
        'results/00_qc/ReadsMultiQCReport.html'
    output:
        'results/00_qc/seqkit.report.raw.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a results/01_cutadapt/*.fastq.gz > {output} '



rule seqkitQCRaw:
    input:
        'results/00_qc/KDRReadsMultiQCReport.html'
    output:
        'results/00_qc/seqkit.report.filtered.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a results/02_kneaddata/*kneaddata.trimmed.fastq > {output} '


