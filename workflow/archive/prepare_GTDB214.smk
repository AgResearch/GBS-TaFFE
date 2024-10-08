# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: 'config/config.yaml'

import os
import pandas as pd


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


wildcard_constraints:
    samples="\w+"

# Global minimum read count for processing
min_reads = 25000
LIBRARY = config["LIBRARY"]

input_fastq_pattern = os.path.join('results', config["LIBRARY"], '01_cutadapt', '{samples}.fastq.gz')
print(input_fastq_pattern)
FIDs, = glob_wildcards(input_fastq_pattern)


rule all:
    input:
        expand('results/{library}/00_QC/seqkit.report.KDR.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.raw.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.bbduk.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.prinseq.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDTrim.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDTRF.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDSILVA138.txt', library = LIBRARY),


checkpoint report_seqkit_raw:
    priority: 1000
    input:
        expand('results/{library}/01_cutadapt/{samples}.fastq.gz', library = LIBRARY, samples = FIDs),
    output:
        os.path.join('results', LIBRARY, '00_QC/seqkit.report.raw.txt')
    benchmark:
        os.path.join('benchmarks', LIBRARY, 'report_seqkit_raw.txt')
    conda:
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 60),
        partition="compute,hugemem"
    shell:
        'seqkit stats -j {threads} -a {input} > {output} '


rule bbduk:
    input:
        reads = 'results/{library}/01_cutadapt/{samples}.fastq.gz',
    output:
        bbdukReads = 'results/{library}/01_readMasking/{samples}.bbduk.fastq.gz'
    log:
        os.path.join('results', '{library}', 'logs', 'bbduk', '{samples}.bbduk.log'),
    benchmark:
        os.path.join('results', '{library}', 'benchmarks', 'bbduk.{samples}.txt'),
    conda:
        'bbduk'
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 8 + ((attempt - 1) * 10),
        partition='compute,hugemem',
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
        'results/{library}/01_readMasking/{samples}.bbduk.fastq.gz'
    output:
        maskedReads = temp('results/{library}/01_readMasking/{samples}.bbduk.prinseq.fastq.gz'),
    log:
        os.path.join('results', '{library}', 'logs', 'prinseq', '{samples}.prinseq.log'),
    benchmark:
        os.path.join('results', '{library}', 'benchmarks', 'prinseq.{samples}.txt'),
    conda:
        'prinseqpp'
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 10),
        partition='compute,hugemem',
    shell:
        'prinseq++ '
        '-threads {threads} '
        '-fastq {input}  '
        '-out_name results/{wildcards.library}/01_readMasking/{wildcards.samples} '
        '-min_len 40 '
        '-lc_entropy=0.5 '
        '-lc_dust=0.5 '
        '-out_gz '
        '2>&1 | tee {log} && '
        'mv results/{wildcards.library}/01_readMasking/{wildcards.samples}_good_out.fastq.gz {output.maskedReads}; '
        'rm results/{wildcards.library}/01_readMasking/{wildcards.samples}_bad_out.fastq.gz '


rule kneaddata:
    input:
        'results/{library}/01_readMasking/{samples}.bbduk.prinseq.fastq.gz'
    output:
        trimReads = temp('results/{library}/02_kneaddata/{samples}.trimmed.fastq'),
        trfReads = temp('results/{library}/02_kneaddata/{samples}.repeats.removed.fastq'),
        silvaReads = temp('results/{library}/02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq'),
        KDRs = temp('results/{library}/02_kneaddata/{samples}.fastq'),
    conda:
        'kneaddata'
    log:
        os.path.join('results', '{library}', 'logs', 'kneaddata', '{samples}.kneaddata.log'),
    benchmark:
        os.path.join('results', '{library}', 'benchmarks', 'kneaddata.{samples}.txt'),
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 20),
        partition='compute,hugemem',
    shell:
        'kneaddata '
        '--trimmomatic-options "ILLUMINACLIP:resources/adapters.fasta:2:30:10 MINLEN:40" '
        '-un {input} '
        '--output-prefix {wildcards.samples} '
        '-t {threads} '
        '--log-level DEBUG '
        '--log {log} '
        '--trimmomatic ~/.conda/envs/kneaddata/share/trimmomatic-0.39-2 ' #TODO Check Path
        '--sequencer-source TruSeq3 '
        '-db /agr/scratch/projects/2023-mbie-rumen-gbs/SILVA138/SILVA_138.1/SLIVA138.1 ' # Embarrassing typo when building index XD
        '-o results/{wildcards.library}/02_kneaddata '
        '&& '
        'touch {output.KDRs} '
        '{output.trimReads} '
        '{output.trfReads} '
        '{output.silvaReads} '


rule gzip_KDR_temps:
    input:
        trimReads = 'results/{library}/02_kneaddata/{samples}.trimmed.fastq',
        trfReads = 'results/{library}/02_kneaddata/{samples}.repeats.removed.fastq',
        silvaReads = 'results/{library}/02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq',
        KDRs ='results/{library}/02_kneaddata/{samples}.fastq',
    output:
        trimReads = temp('results/{library}/02_kneaddata/{samples}.trimmed.fastq.gz'),
        trfReads = temp('results/{library}/02_kneaddata/{samples}.repeats.removed.fastq.gz'),
        silvaReads = temp('results/{library}/02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq.gz'),
        KDRs ='results/{library}/02_kneaddata/{samples}.fastq.gz',
    benchmark:
        os.path.join('results', '{library}', 'benchmarks', 'gzip_KDR_temps.{samples}.txt'),
    conda:
        "pigz"
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 16 + ((attempt - 1) * 16),
        partition='compute,hugemem',
    shell:
        """
        
        pigz -p {threads} -c {input.KDRs} > {output.KDRs} && 
        pigz -p {threads} -c {input.trimReads} > {output.trimReads} && 
        pigz -p {threads} -c {input.trfReads} > {output.trfReads} && 
        pigz -p {threads} -c {input.silvaReads} > {output.silvaReads} &&

        rm {input.KDRs} {input.trimReads} {input.trfReads} {input.silvaReads};

        """


def get_seqkitKneaddata_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "02_kneaddata/{samples}.fastq.gz"), samples = passed)


rule report_seqkit_KDR:
    priority: 1000
    input:
        KDRs = get_seqkitKneaddata_passing_samples,
    output:
        report = os.path.join('results', LIBRARY, '00_QC' ,'seqkit.report.KDR.txt'),
        token = os.path.join('results', LIBRARY, '00_QC' ,'.priority.semaphore.tkn')
    benchmark:
        os.path.join('results', LIBRARY, 'benchmarks', 'report_seqkit_KDR.txt'),
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute,hugemem',
    shell:
        'seqkit stats -j {threads} -a {input.KDRs} > {output.report} && '
        'touch {output.token} '


def get_seqkitMaskingBBDukReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "01_readMasking/{samples}.bbduk.fastq.gz"), samples = passed)


rule report_seqkit_bbduk:
    input:
        bbdukReads = get_seqkitMaskingBBDukReads_passing_samples,
        token = os.path.join('results', LIBRARY, '00_QC' ,'.priority.semaphore.tkn'),
    output:
        os.path.join("results", LIBRARY, "00_QC", "seqkit.report.bbduk.txt")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "report_seqkit_bbduk.txt"),
    conda:
        'envs/seqkit.yaml'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 90 + ((attempt - 1) * 30),
        partition='compute,hugemem'
    shell:
        'seqkit stats -j {threads} -a {input.bbdukReads} > {output} '


def get_seqkitMaskingPrinseqReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "01_readMasking/{samples}.bbduk.prinseq.fastq.gz"), samples = passed)


rule report_seqkit_prinseq:
    input:
        prinseqReads = get_seqkitMaskingPrinseqReads_passing_samples,
        token = os.path.join('results', LIBRARY, '00_QC' ,'.priority.semaphore.tkn'),
    output:
        os.path.join("results", LIBRARY, "00_QC", "seqkit.report.prinseq.txt")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "report_seqkit_prinseq.txt"),
    conda:
        'envs/seqkit.yaml'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute,hugemem',
    shell:
        'seqkit stats -j {threads} -a {input.prinseqReads} > {output} '


def get_seqkitKneaddataTrimReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "02_kneaddata/{samples}.trimmed.fastq.gz"), samples = passed)


rule report_seqkit_KDTrim: 
    input:
        trimReads = get_seqkitKneaddataTrimReads_passing_samples,
        token = os.path.join('results', LIBRARY, '00_QC' ,'.priority.semaphore.tkn'),
    output:
        os.path.join('results', LIBRARY, '00_QC', 'seqkit.report.KDTrim.txt')
    benchmark:
        os.path.join('results', LIBRARY, 'benchmarks', 'report_seqkit_KDSILVA.txt'),
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute,hugemem',
    shell:
        'seqkit stats -j {threads} -a {input.trimReads} > {output} '


def get_seqkitKneaddataTRFReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "02_kneaddata/{samples}.repeats.removed.fastq.gz"), samples = passed)


rule report_seqkit_KDTRF:
    input:
        trfReads = get_seqkitKneaddataTRFReads_passing_samples,
        token = os.path.join('results', LIBRARY, '00_QC' ,'.priority.semaphore.tkn'),
    output:
        os.path.join('results', LIBRARY, '00_QC', 'seqkit.report.KDTRF.txt')
    benchmark:
        os.path.join('results', LIBRARY, 'benchmarks', 'report_seqkit_KDSILVA.txt')
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute,hugemem',
    shell:
        'seqkit stats -j {threads} -a {input.trfReads} > {output} '


def get_seqkitKneaddataSILVAReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq.gz"), samples = passed)


rule report_seqkit_KDSILVA138:
    input:
        silvaReads = get_seqkitKneaddataSILVAReads_passing_samples,
        token = os.path.join('results', LIBRARY, '00_QC' ,'.priority.semaphore.tkn'),
    output:
        os.path.join('results', LIBRARY, '00_QC', 'seqkit.report.KDSILVA138.txt')
    benchmark:
        os.path.join('results', LIBRARY, 'benchmarks', 'report_seqkit_KDSILVA.txt')
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute,hugemem',
    shell:
        'seqkit stats -j {threads} -a {input.silvaReads} > {output} '
