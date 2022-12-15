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
        expand('03_kraken2/{samples}.report.k2', samples = FIDs),
        expand('03_kraken2GTDB/{samples}.GTDB.report.k2', samples = FIDs),
        # expand('03_metaphlan4/{samples}.metaphlan4.profile.txt', samples = FIDs),
        # expand('03_humann/{samples}_kneaddata_pathabundance.tsv', samples = FIDs),
        '00_qc/ReadsMultiQCReport.html',
        '00_qc/KDRReadsMultiQCReport.html'


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
        'Spewing barcodes for: {params.libraries}...\n'
    shell:
        'gquery '
        '-t gbs_keyfile '
        '-b library -p "columns=factid,barcode;fasta;noheading;no_unpivot" '
        '{params.libraries} > '
        '{output.barcodes} 2> '
        '{log}'


rule cutadapt: #demultiplexing GBS reads
    output:
        expand('01_cutadapt/{samples}.fastq.gz', samples = FIDs),
    input:
        barcodes = rules.generateBarcodes.output.barcodes,
        lane01 = config['novaseq']['lane01'],
        lane02 = config['novaseq']['lane02'],
    conda:
        'cutadapt'
        # 'docker://quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'
    threads: 16
    resources:
        mem_gb=32,
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
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
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
    conda:
        'multiqc'
        # 'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '-n 00_qc/ReadsMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'



#TODO Rule to Build Rambv2 index



rule kneaddata:
    output:
        trimReads = temp('02_kneaddata/{samples}_kneaddata.trimmed.fastq'),
        trfReads = temp('02_kneaddata/{samples}_kneaddata.repeats.removed.fastq'),
        ovineReads = temp('02_kneaddata/{samples}_kneaddata_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq'),
        clnReads = '02_kneaddata/{samples}_kneaddata.fastq',
        readStats = '02_kneaddata/{samples}.read.stats.txt'
    input:
        reads = '01_cutadapt/{samples}.fastq.gz',
    conda:
        'biobakery'
    log:
        'logs/{samples}.kneaddata.log'
    threads: 8
    resources:
        mem_gb=12,
        time='02:00:00'
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
        "00_qc/fastqcKDR/{samples}_kneaddata_fastqc.zip"
    input:
        fastq = '02_kneaddata/{samples}_kneaddata.fastq'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
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
        '00_qc/KDRReadsMultiQCReport.html'
    input:
        fastqc= expand('00_qc/fastqcKDR/{samples}_kneaddata_fastqc.zip', samples = FIDs)
    conda:
        'multiqc'
        # 'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '-n 00_qc/KDRReadsMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'



# TODO derepliation ???



rule metaphlan4:
    output:
        '03_metaphlan4/{samples}.metaphlan4.profile.txt'
    input:
        KDRs=rules.kneaddata.output.clnReads,
    conda:
        'metaphlan'
    log:
        'logs/{samples}.metaphlan4.log'
    threads: 4
    message:
        'Running Metaphlan4 on: {wildcards.samples} \n'
    shell:
        'metaphlan '
        '{input} '
        '--input_type fastq '
        '--bowtie2db ref/biobakery/metaphlan '
        '--nproc {threads} '
        #'--nreads $(cat 02_kneaddata/{input.KDRs} | grep "^+$" | wc -l) ' #TODO update to zcat when using compressed reads
        '--unclassified_estimation '
        '-t rel_ab '
        '> {output} '



rule kraken2:
    output:
        k2Out='03_kraken2/{samples}.out.k2',
        k2Report='03_kraken2/{samples}.report.k2'
    input:
        KDRs=rules.kneaddata.output.clnReads
    log:
        'logs/{samples}.kraken2.log'
    conda:
        'kraken2'
    threads: 12 
    resources: 
        mem_gb=180,
        partition="inv-bigmem,inv-bigmem-fast"
    shell:
        'kraken2 '
        '--db ref/kraken2 '
        '--report {output.k2Report} '
        '--report-minimizer-data '
        '{input.KDRs} > {output.k2Out}'



rule kraken2GTDB:
    output:
        k2OutGTDB='03_kraken2GTDB/{samples}.GTDB.out.k2',
        k2ReportGTDB='03_kraken2GTDB/{samples}.GTDB.report.k2'
    input:
        KDRs=rules.kneaddata.output.clnReads
    log:
        'logs/{samples}.kraken2.GTDB.log'
    conda:
        'kraken2'
    threads: 18 
    resources: 
        mem_gb=330,
        partition="inv-bigmem,inv-bigmem-fast"
    shell:
        'kraken2 '
        '--db /dataset/2022-BJP-GTDB/active/kraken/GTDB '
        '--report {output.k2ReportGTDB} '
        '--report-minimizer-data '
        '{input.KDRs} > {output.k2OutGTDB}'



rule braken:
    output:
        '03_kraken2/{samples}.braken.out'
    input:
        k2out=rules.kraken2.output
    log:
        'logs/{samples}.braken.log'
    conda:
        ''# TBD
    threads: 8
    message:
        ''
    shell:
        ''



rule humann3:
    output:
        genes = '03_humann/{samples}_kneaddata_genefamilies.tsv',
        pathways = '03_humann/{samples}_kneaddata_pathabundance.tsv',
        pathwaysCoverage = '03_humann/{samples}_kneaddata_pathcoverage.tsv'
    input:
        KDRs=rules.kneaddata.output.clnReads
    log:
        'logs/{samples}.human3.log'
    conda:
        'biobakery'
    threads: 8
    resources:
        mem_gb=24,
    message:
        'humann3 profiling: {wildcards.samples}\n'
    shell:
        'humann '
        '--threads {threads} '
        '--input {input.KDRs} '
        '--output 03_humann '
        '--nucleotide-database ref/biobakery/humann_dbs/chocophlan '
        '--bypass-prescreen '
        #'--bypass-nucleotide-search '
        '--memory-use maximum '
        '--input-format fastq '
        '--search-mode uniref90 '
        '--protein-database ref/biobakery/humann_dbs/uniref ' #TODO Update
        '--verbose '
        '--log-level INFO '
        '--o-log {log} '
        '--remove-temp-output '



# # rule human3GTDB:
#     output:

#     input:

#     log:

#     threads:

#     message:

#     shell:

