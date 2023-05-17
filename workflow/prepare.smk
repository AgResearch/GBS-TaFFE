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

	
rule all:
    input:
        'results/00_QC/ReadsMultiQCReport.html',
        'results/00_QC/seqkit.report.KDTrim.txt',
        'results/00_QC/seqkit.report.KDTRF.txt',
        'results/00_QC/seqkit.report.KDOvine.txt',
        'results/00_QC/seqkit.report.KDSILVA138.txt',
        'results/00_QC/seqkit.report.raw.txt',
        'results/00_QC/seqkit.report.bbduk.txt',
        'results/00_QC/seqkit.report.prinseq.txt',
        'results/00_QC/seqkit.report.KDR.txt',


# RULES TO GENERATE BARCODES AND DEMULTIPLEX
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
        run = config['novaseq']
    output:
        expand('results/01_cutadapt/{samples}.fastq.gz', samples = FIDs),
    conda:
        'cutadapt'
        # 'docker://quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'
    threads: 32
    resources:
        mem_gb=64,
        time="12:00:00"
    message:
        'Demultiplexing lanes...'
    shell:
        'zcat {input.run} | '
        'cutadapt '
        '-j {threads} '
        '--discard-untrimmed '
        '--no-indels '
        '-g ^file:{input.barcodes} '
        r'-o "results/01_cutadapt/{{name}}.fastq.gz" '
        '-' # indicates stdin

# STANDARD READ FILTERING AND QC RULES
rule fastqc:
    input:
        fastq = 'downsampled/{sample}.fastq.gz'
    output:
        html = 'results/00_QC/fastqc/{sample}_fastqc.html',
        zip = 'results/00_QC/fastqc/{sample}_fastqc.zip'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
    message:
        'Running QC on reads: {wildcards.sample}\n'
    shell:
        'fastqc '
        '-o results/00_QC/fastqc/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule bbduk:
    input:
        reads = 'downsampled/{sample}.fastq.gz',
    output:
        bbdukReads = 'results/01_readMasking/{sample}.bbduk.fastq.gz'
    log:
        'logs/bbduk/{sample}.bbduk.log'
    conda:
        'bbduk'
    threads:4
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.reads} '
        'entropy=0.3 '
        'entropywindow=50 '
        'trimpolygright=5 '
        'qtrim=r '
        'trimq=20 '
        'out={output.bbdukReads} '
        '2>&1 | tee {log}'


rule prinseq:
    input:
        'results/01_readMasking/{sample}.bbduk.fastq.gz'
    output:
        maskedReads = 'results/01_readMasking/{sample}.bbduk.prinseq.fastq.gz',
        badReads = temp('results/01_readMasking/{sample}_bad_out.fastq.gz'),
    log:
        'logs/prinseq/{sample}.prinseq.log'
    conda:
        'prinseqPP'
    threads:4
    shell:
        'prinseq++ '
        '-threads {threads} '
        '-fastq {input}  '
        '-out_name results/01_readMasking/{wildcards.sample} '
        '-min_len 40 '
        '-lc_entropy=0.5 '
        '-lc_dust=0.5 '
        '-out_gz '
        '2>&1 | tee {log} && '
        'mv results/01_readMasking/{wildcards.sample}_good_out.fastq.gz {output.maskedReads} '


rule kneaddata:
    input:
        'results/01_readMasking/{sample}.bbduk.prinseq.fastq.gz'
    output:
        trimReads = temp('results/02_kneaddata/{sample}.trimmed.fastq'),
        trfReads = temp('results/02_kneaddata/{sample}.repeats.removed.fastq'),
        ovineReads = temp('results/02_kneaddata/{sample}_ARS_UCD1.3_bowtie2_contam.fastq'),
        silvaReads = temp('results/02_kneaddata/{sample}_SLIVA138.1_bowtie2_contam.fastq'),
        KDRs ='results/02_kneaddata/{sample}.fastq',
    conda:
        'biobakery'
    log:
        'logs/kneaddata/{sample}.kneaddata.log'
    threads: 8
    resources:
        mem_gb=8,
        time='02:00:00'
    message:
        'kneaddata: {wildcards.sample}\n'
    shell:
        'kneaddata '
        '--trimmomatic-options "ILLUMINACLIP:/home/perrybe/conda-envs/biobakery/share/trimmomatic-0.39-2/adapters/illuminaAdapters.fa:2:30:10 MINLEN:40" '
        '--input {input} '
        '--output-prefix {wildcards.sample} '
        '-t {threads} '
        '--log-level DEBUG '
        #'--log {log} '
        '--trimmomatic /home/perrybe/conda-envs/biobakery/share/trimmomatic '
        '--sequencer-source TruSeq3 '
        '-db /home/perrybe/quickQC/cow_rumens/SMK-PROFILE/ref/ARS_UCD1.3 '
        '-db /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/SILVA_138.1/SLIVA138.1 ' # Embarrassing typo when building index XD
        '-o results/02_kneaddata '
        '2>&1 | tee {log} '


rule fastqcKDRs:
    input:
        fastq = 'results/02_kneaddata/{sample}.fastq'
    output:
        'results/00_QC/fastqcKDR/{sample}_fastqc.zip'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
    message:
        'Running QC on reads: {wildcards.sample}\n'
    shell:
        'fastqc '
        '-o results/00_QC/fastqcKDR/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule seqkitKneaddataTrimReads: #TODO expand these
    input:
        trimReads = expand('results/02_kneaddata/{sample}.trimmed.fastq', sample = FIDs),
    output:
        'results/00_QC/seqkit.report.KDTrim.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a {input.trimReads} > {output} '


rule seqkitKneaddataTRFReads:
    input:
        trfReads = expand('results/02_kneaddata/{sample}.repeats.removed.fastq', sample = FIDs),
    output:
        'results/00_QC/seqkit.report.KDTRF.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a {input.trfReads} > {output} '


rule seqkitKneaddataHostReads:
    input:
        HostReads = expand('results/02_kneaddata/{sample}_ARS_UCD1.3_bowtie2_contam.fastq', sample = FIDs),
    output:
        'results/00_QC/seqkit.report.KDOvine.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a {input.ovineReads} > {output} '


rule seqkitKneaddataSILVAReads:
    input:
        silvaReads = expand('results/02_kneaddata/{sample}_SLIVA138.1_bowtie2_contam.fastq', sample = FIDs),
    output:
        'results/00_QC/seqkit.report.KDSILVA138.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a {input.silvaReads} > {output} '


rule seqkitRaw:
    output:
        'results/00_QC/seqkit.report.raw.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a fastq/*.fastq.gz > {output} '


rule seqkitMaskingBBDukReads:
    input:
        bbdukReads = expand('results/01_readMasking/{sample}.bbduk.fastq.gz', sample = FIDs),
    output:
        'results/00_QC/seqkit.report.bbduk.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a {input.bbdukReads} > {output} '


rule seqkitMaskingPrinseqReads:
    input:
        prinseqReads = expand('results/01_readMasking/{sample}.bbduk.prinseq.fastq.gz', sample = FIDs),
    output:
        'results/00_QC/seqkit.report.prinseq.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a {input.prinseqReads} > {output} '


rule seqkitKneaddata:
    input:
        KDRs = expand('results/02_kneaddata/{sample}.fastq', sample = FIDs),
    output:
        'results/00_QC/seqkit.report.KDR.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a {input.KDRs} > {output} '