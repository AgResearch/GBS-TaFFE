# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: 'config/config.yaml'

import os
import pandas as pd


FIDs, = glob_wildcards('results/01_cutadapt/{samples}.fastq.gz')


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
        'results/00_QC/seqkit.report.KDTrim.txt',
        'results/00_QC/seqkit.report.KDTRF.txt',
        'results/00_QC/seqkit.report.KDOvis.txt',
        'results/00_QC/seqkit.report.KDBos.txt',
        'results/00_QC/seqkit.report.KDCapra.txt',
        'results/00_QC/seqkit.report.KDCervus.txt',
        'results/00_QC/seqkit.report.KDSILVA138.txt',
        'results/00_QC/seqkit.report.raw.txt',
        'results/00_QC/seqkit.report.bbduk.txt',
        'results/00_QC/seqkit.report.prinseq.txt',
        'results/00_QC/seqkit.report.KDR.txt',


checkpoint seqkitRaw:
    input:
        expand('results/01_cutadapt/{samples}.fastq.gz', samples = FIDs),
    output:
        'results/00_QC/seqkit.report.raw.txt'
    benchmark:
        'benchmarks/seqkitRaw.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0' 
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=2,
        time="00:60:00",
        partition="large,milan"
    threads: 12
    shell:
        'seqkit stats -j {threads} -a {input} > {output} '


# STANDARD READ FILTERING AND QC RULES
rule bbduk:
    input:
        reads = 'results/01_cutadapt/{samples}.fastq.gz',
    output:
        bbdukReads = temp('results/01_readMasking/{samples}.bbduk.fastq.gz')
    log:
        'logs/bbduk/{samples}.bbduk.log'
#    benchmark:
#        'benchmarks/bbduk.{samples}.txt'
    #container:
    #    'docker://quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' 
    conda:
        #'env/bbduk.yaml'
        'bbduk'
    threads:4
    resources:
        mem_gb=1,
        time="00:01:00",
        partition='large,milan',
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
        'results/01_readMasking/{samples}.bbduk.fastq.gz'
    output:
        maskedReads = temp('results/01_readMasking/{samples}.bbduk.prinseq.fastq.gz'),
        badReads = temp('results/01_readMasking/{samples}_bad_out.fastq.gz'),
    log:
        'logs/prinseq/{samples}.prinseq.log'
#    benchmark:
#        'benchmarks/prinseq.{samples}.txt'
    #container:
    #    'docker://quay.io/biocontainers/prinseq-plus-plus:1.2.4--h7ff8a90_2' 
    conda:
        #'env/prinseqPP.yaml'
        'prinseqpp'
    resources:
        mem_gb=1,
        time="00:04:00",
        partition='large,milan',
    threads:4
    shell:
        'prinseq++ '
        '-threads {threads} '
        '-fastq {input}  '
        '-out_name results/01_readMasking/{wildcards.samples} '
        '-min_len 40 '
        '-lc_entropy=0.5 '
        '-lc_dust=0.5 '
        '-out_gz '
        '2>&1 | tee {log} && '
        'mv results/01_readMasking/{wildcards.samples}_good_out.fastq.gz {output.maskedReads} '


rule sana:
    input:
        maskedReads = 'results/01_readMasking/{samples}.bbduk.prinseq.fastq.gz'
    output:
        sanitised = temp('results/01_readMasking/{samples}.bbduk.prinseq.sana.fastq.gz'),
    log:
        'logs/sasna.{samples}.log'
#    benchmark:
#        'benchmarks/prinseq.{samples}.txt'
    #container:
    #    'docker://quay.io/biocontainers/prinseq-plus-plus:1.2.4--h7ff8a90_2' 
    conda:
        #'env/prinseqPP.yaml'
        'seqkit'
    resources:
        mem_gb=1,
        time="00:01:00",
        partition='large,milan',
    threads:4
    shell:
        'seqkit sana {input.maskedReads} > {output.sanitised} '
        '2> {log}'


def get_seqkitMaskingPrinseqReads_passing_samples(wildcards):
    file = checkpoints.seqkitRaw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/01_readMasking/{samples}.bbduk.prinseq.sana.fastq.gz", samples = passed)

checkpoint seqkitMaskingPrinseqReads:
    input:
        prinseqReads = get_seqkitMaskingPrinseqReads_passing_samples,
    output:
        'results/00_QC/seqkit.report.prinseq.txt'
    benchmark:
        'benchmarks/seqkitMaskingPrinseqReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.prinseqReads} > {output} '


rule kneaddata:
    input:
        'results/01_readMasking/{samples}.bbduk.prinseq.fastq.gz'
    output:
        trimReads = temp('results/02_kneaddata/{samples}.trimmed.fastq'),
        trfReads = temp('results/02_kneaddata/{samples}.repeats.removed.fastq'),
        ovineReads = temp('results/02_kneaddata/{samples}_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq'),
        bosReads = temp('results/02_kneaddata/{samples}_ARS_UCD1.3_bowtie2_contam.fastq'),
        capraReads = temp('results/02_kneaddata/{samples}_CAPRA_ARS1.2_bowtie2_contam.fastq'),
        cervusReads = temp('results/02_kneaddata/{samples}_mCerEla1_bowtie2_contam.fastq'),
        silvaReads = temp('results/02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq'),
        KDRs ='results/02_kneaddata/{samples}.fastq',
    #container:
    #    'docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0'
    conda:
        #'env/biobakery.yaml'
        'kneaddata'
    log:
        'logs/kneaddata/{samples}.kneaddata.log'
#    benchmark:
#        'benchmarks/kneaddata.{samples}.txt'
    threads: 8
    resources:
        mem_gb=24,
        time="00:14:00",
        partition='large,milan',
    message:
        'kneaddata: {wildcards.samples}\n'
    shell:
        'kneaddata '
        '--trimmomatic-options "ILLUMINACLIP:/nesi/nobackup/agresearch03843/illuminaAdapters.fa:2:30:10 MINLEN:40" '
        '-un {input} '
        '--output-prefix {wildcards.samples} '
        '-t {threads} '
        '--log-level DEBUG '
        '--log {log} '
        '--trimmomatic /home/perrybe/.conda/envs/kneaddata/share/trimmomatic-0.39-2 '
        '--sequencer-source TruSeq3 '
        '-db /nesi/nobackup/agresearch03843/RAMBV2/Rambv2/GCF_016772045.1-ARS-UI-Ramb-v2.0 '
        '-db /nesi/nobackup/agresearch03843/ARSUCD1/ARS_UCD1.3 '
        '-db /nesi/nobackup/agresearch03843/CAPRA/CAPRA_ARS1.2 '
        '-db /nesi/nobackup/agresearch03843/CERVUS/mCerEla1 '
        '-db /nesi/nobackup/agresearch03843/SILVA138/SILVA_138.1/SLIVA138.1 ' # Embarrassing typo when building index XD
        '-o results/02_kneaddata '


def get_seqkitKneaddataTrimReads_passing_samples(wildcards):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/02_kneaddata/{samples}.trimmed.fastq", samples = passed)


rule seqkitKneaddataTrimReads: 
    input:
        trimReads = get_seqkitKneaddataTrimReads_passing_samples
    output:
        'results/00_QC/seqkit.report.KDTrim.txt'
    benchmark:
        'benchmarks/seqkitKneaddataTrimReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.trimReads} > {output} '


def get_seqkitKneaddataTRFReads_passing_samples(wildcards):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/02_kneaddata/{samples}.repeats.removed.fastq", samples = passed)


rule seqkitKneaddataTRFReads:
    input:
        trfReads = get_seqkitKneaddataTRFReads_passing_samples
    output:
        'results/00_QC/seqkit.report.KDTRF.txt'
    benchmark:
        'benchmarks/seqkitKneaddataTRFReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.trfReads} > {output} '


def get_seqkitKneaddataOvisReads_passing_samples(wildcards):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/02_kneaddata/{samples}_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq", samples = passed)


rule seqkitKneaddataOvisReads:
    input:
        HostReads = get_seqkitKneaddataOvisReads_passing_samples,
    output:
        'results/00_QC/seqkit.report.KDOvis.txt'
    benchmark:
        'benchmarks/seqkitKneaddataHostReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.HostReads} > {output} '

def get_seqkitKneaddataBosReads_passing_samples(wildcards):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/02_kneaddata/{samples}_ARS_UCD1.3_bowtie2_contam.fastq", samples = passed)


rule seqkitKneaddataBosReads:
    input:
        HostReads = get_seqkitKneaddataBosReads_passing_samples,
    output:
        'results/00_QC/seqkit.report.KDBos.txt'
    benchmark:
        'benchmarks/seqkitKneaddataHostReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.HostReads} > {output} '


def get_seqkitKneaddataCapraReads_passing_samples(wildcards):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/02_kneaddata/{samples}_CAPRA_ARS1.2_bowtie2_contam.fastq", samples = passed)


rule seqkitKneaddataCapraReads:
    input:
        HostReads = get_seqkitKneaddataCapraReads_passing_samples,
    output:
        'results/00_QC/seqkit.report.KDCapra.txt'
    benchmark:
        'benchmarks/seqkitKneaddataCapraReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.HostReads} > {output} '

def get_seqkitKneaddataCervusReads_passing_samples(wildcards):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/02_kneaddata/{samples}_mCerEla1_bowtie2_contam.fastq", samples = passed)


rule seqkitKneaddataCervusReads:
    input:
        HostReads = get_seqkitKneaddataCervusReads_passing_samples,
    output:
        'results/00_QC/seqkit.report.KDCervus.txt'
    benchmark:
        'benchmarks/seqkitKneaddataCervusReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.HostReads} > {output} '


def get_seqkitKneaddataSILVAReads_passing_samples(wildcards):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq", samples = passed)


rule seqkitKneaddataSILVAReads:
    input:
        silvaReads = get_seqkitKneaddataSILVAReads_passing_samples,
    output:
        'results/00_QC/seqkit.report.KDSILVA138.txt'
    benchmark:
        'benchmarks/seqkitKneaddataSILVAReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.silvaReads} > {output} '


def get_seqkitMaskingBBDukReads_passing_samples(wildcards):
    file = checkpoints.seqkitRaw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/01_readMasking/{samples}.bbduk.fastq.gz", samples = passed)


rule seqkitMaskingBBDukReads:
    input:
        bbdukReads = get_seqkitMaskingBBDukReads_passing_samples,
    output:
        'results/00_QC/seqkit.report.bbduk.txt'
    benchmark:
        'benchmarks/seqkitMaskingBBDukReads.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.bbdukReads} > {output} '


def get_seqkitKneaddata_passing_samples(wildcards):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/02_kneaddata/{samples}.fastq", samples = passed)


rule seqkitKneaddata:
    input:
        KDRs = get_seqkitKneaddata_passing_samples,
    output:
        'results/00_QC/seqkit.report.KDR.txt'
    benchmark:
        'benchmarks/seqkitKneaddata.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    resources:
        mem_gb=4,
        time="00:90:00",
        partition="large,milan"
    threads: 32
    shell:
        'seqkit stats -j {threads} -a {input.KDRs} > {output} '
