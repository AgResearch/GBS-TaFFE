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
LIBRARY, = config["LIBRARY"]

input_fastq_pattern = os.path.join('results', config["LIBRARY"], '01_cutadapt', '{samples}.fastq.gz')
print(input_fastq_pattern)
FIDs, = glob_wildcards(input_fastq_pattern)


rule all:
    input:
        expand('results/{library}/00_QC/seqkit.report.raw.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.prinseq.txt', library = LIBRARY),
        expand("results/{library}/kraken2.{library}.genus.counts.tsv",  library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDTrim.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDTRF.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDOvis.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDBos.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDCapra.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDCervus.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDSILVA138.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.bbduk.txt', library = LIBRARY),


checkpoint seqkitRaw:
    priority: 1000
    input:
        expand('results/{library}/01_cutadapt/{samples}.fastq.gz', library = LIBRARY, samples = FIDs),
    output:
        os.path.join('results', LIBRARY, '00_QC/seqkit.report.raw.txt')
    benchmark:
        os.path.join('benchmarks', LIBRARY, 'seqkitRaw.txt')
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 60),
        partition="compute"
    shell:
        'seqkit stats -j {threads} -a {input} > {output} '


rule bbduk:
    input:
        reads = 'results/{library}/01_cutadapt/{samples}.fastq.gz',
    output:
        bbdukReads = temp('results/{library}/01_readMasking/{samples}.bbduk.fastq.gz')
    log:
        'logs/{library}/bbduk/{samples}.bbduk.log'
    conda:
        'bbduk'
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 8 + ((attempt - 1) * 10),
        partition='compute',
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
        badReads = temp('results/{library}/01_readMasking/{samples}_bad_out.fastq.gz'),
    log:
        'logs/{library}/prinseq/{samples}.prinseq.log'
    conda:
        'prinseqpp'
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 10 + ((attempt - 1) * 10),
        partition='compute',
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
        'mv results/{wildcards.library}/01_readMasking/{wildcards.samples}_good_out.fastq.gz {output.maskedReads} '


def get_seqkitMaskingPrinseqReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitRaw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "/01_readMasking/{samples}.bbduk.prinseq.fastq.gz"), samples = passed)


checkpoint seqkitMaskingPrinseqReads:
    input:
        prinseqReads = get_seqkitMaskingPrinseqReads_passing_samples,
    output:
        os.path.join("results", LIBRARY, "00_QC", "seqkit.report.prinseq.txt")
    benchmark:
        os.path.join("benchmarks", LIBRARY, "seqkitMaskingPrinseqReads.txt")
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.prinseqReads} > {output} '


def get_seqkitMaskingBBDukReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitRaw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/01_readMasking/{samples}.bbduk.fastq.gz", samples = passed, library=lib)


rule seqkitMaskingBBDukReads:
    input:
        bbdukReads = get_seqkitMaskingBBDukReads_passing_samples,
    output:
        'results/{library}/00_QC/seqkit.report.bbduk.txt'
    benchmark:
        'benchmarks/{library}/seqkitMaskingBBDukReads.txt'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 90 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.bbdukReads} > {output} '


rule kneaddata:
    input:
        'results/{library}/01_readMasking/{samples}.bbduk.prinseq.fastq.gz'
    output:
        trimReads = temp('results/{library}/02_kneaddata/{samples}.trimmed.fastq'),
        trfReads = temp('results/{library}/02_kneaddata/{samples}.repeats.removed.fastq'),
        ovineReads = temp('results/{library}/02_kneaddata/{samples}_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq'),
        bosReads = temp('results/{library}/02_kneaddata/{samples}_ARS_UCD1.3_bowtie2_contam.fastq'),
        capraReads = temp('results/{library}/02_kneaddata/{samples}_CAPRA_ARS1.2_bowtie2_contam.fastq'),
        cervusReads = temp('results/{library}/02_kneaddata/{samples}_mCerEla1_bowtie2_contam.fastq'),
        silvaReads = temp('results/{library}/02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq'),
        KDRs = temp('results/{library}/02_kneaddata/{samples}.fastq'),
    conda:
        'kneaddata'
    log:
        'logs/{library}/kneaddata/{samples}.kneaddata.log'
    benchmark:
        'benchmarks/{library}/kneaddata.{samples}.txt'
    threads: 18
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 20),
        partition='compute',
    shell:
        'kneaddata '
        '--trimmomatic-options "ILLUMINACLIP:resources/illuminaAdapters.fa:2:30:10 MINLEN:40" '
        '-un {input} '
        '--output-prefix {wildcards.samples} '
        '-t {threads} '
        '--log-level DEBUG '
        '--log {log} '
        '--trimmomatic ~/.conda/envs/kneaddata/share/trimmomatic-0.39-2 '
        '--sequencer-source TruSeq3 '
        '-db /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/Rambv2/GCF_016772045.1-ARS-UI-Ramb-v2.0 '
        '-db /nesi/nobackup/agresearch03843/ARSUCD1/ARS_UCD1.3 '
        '-db /nesi/nobackup/agresearch03843/CAPRA/CAPRA_ARS1.2 '
        '-db /nesi/nobackup/agresearch03843/CERVUS/mCerEla1 '
        '-db /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/SILVA_138.1/SLIVA138.1 ' # Embarrassing typo when building index XD
        '-o results/{wildcards.library}/02_kneaddata '
        '&& '
        'touch {output.KDRs} '
        '{output.trimReads} '
        '{output.trfReads} '
        '{output.ovineReads} '
        '{output.bosReads} '
        '{output.capraReads} '
        '{output.cervusReads} '
        '{output.silvaReads} '


rule gzip_KDR_temps:
    input:
        trimReads = 'results/{library}/02_kneaddata/{samples}.trimmed.fastq',
        trfReads = 'results/{library}/02_kneaddata/{samples}.repeats.removed.fastq',
        ovineReads = 'results/{library}/02_kneaddata/{samples}_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq',
        bosReads = 'results/{library}/02_kneaddata/{samples}_ARS_UCD1.3_bowtie2_contam.fastq',
        capraReads = 'results/{library}/02_kneaddata/{samples}_CAPRA_ARS1.2_bowtie2_contam.fastq',
        cervusReads = 'results/{library}/02_kneaddata/{samples}_mCerEla1_bowtie2_contam.fastq',
        silvaReads = 'results/{library}/02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq',
        KDRs ='results/{library}/02_kneaddata/{samples}.fastq',
    output:
        trimReads = temp('results/{library}/02_kneaddata/{samples}.trimmed.fastq.gz'),
        trfReads = temp('results/{library}/02_kneaddata/{samples}.repeats.removed.fastq.gz'),
        ovineReads = 'results/{library}/02_kneaddata/{samples}_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq.gz',
        bosReads = temp('results/{library}/02_kneaddata/{samples}_ARS_UCD1.3_bowtie2_contam.fastq.gz'),
        capraReads = temp('results/{library}/02_kneaddata/{samples}_CAPRA_ARS1.2_bowtie2_contam.fastq.gz'),
        cervusReads = temp('results/{library}/02_kneaddata/{samples}_mCerEla1_bowtie2_contam.fastq.gz'),
        silvaReads = temp('results/{library}/02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq.gz'),
        KDRs ='results/{library}/02_kneaddata/{samples}.fastq.gz',
    conda:
        "pigz"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 15 + ((attempt - 1) * 12),
        partition='compute',
    shell:
        "pigz -p 12 "
        "{input.KDRs} "
        "{input.trimReads} "
        "{input.trfReads} "
        "{input.ovineReads} "
        "{input.bosReads} "
        "{input.capraReads} "
        "{input.cervusReads} "
        "{input.silvaReads} "


def get_seqkitKneaddata_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/02_kneaddata/{samples}.fastq.gz", samples = passed, library = lib)


rule seqkitKneaddata:
    input:
        KDRs = get_seqkitKneaddata_passing_samples,
    output:
        'results/{library}/00_QC/seqkit.report.KDR.txt'
    benchmark:
        'benchmarks/{library}/seqkitKneaddata.txt'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.KDRs} > {output} '


def get_seqkitKneaddataTrimReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/02_kneaddata/{samples}.trimmed.fastq.gz", samples = passed, library = lib)


rule seqkitKneaddataTrimReads: 
    input:
        trimReads = get_seqkitKneaddataTrimReads_passing_samples
    output:
        'results/{library}/00_QC/seqkit.report.KDTrim.txt'
    benchmark:
        'benchmarks/{library}/seqkitKneaddataTrimReads.txt'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.trimReads} > {output} '


def get_seqkitKneaddataTRFReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/02_kneaddata/{samples}.repeats.removed.fastq.gz", samples = passed, library = lib)


rule seqkitKneaddataTRFReads:
    input:
        trfReads = get_seqkitKneaddataTRFReads_passing_samples
    output:
        'results/{library}/00_QC/seqkit.report.KDTRF.txt'
    benchmark:
        'benchmarks/{library}/seqkitKneaddataTRFReads.txt'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.trfReads} > {output} '


def get_seqkitKneaddataOvisReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/02_kneaddata/{samples}_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq.gz", samples = passed, library = lib)


rule seqkitKneaddataOvisReads:
    input:
        HostReads = get_seqkitKneaddataOvisReads_passing_samples,
    output:
        'results/{library}/00_QC/seqkit.report.KDOvis.txt'
    benchmark:
        'benchmarks/{library}/seqkitKneaddataHostReads.txt'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.HostReads} > {output} '

def get_seqkitKneaddataBosReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/02_kneaddata/{samples}_ARS_UCD1.3_bowtie2_contam.fastq.gz", samples = passed, library = lib)


rule seqkitKneaddataBosReads:
    input:
        HostReads = get_seqkitKneaddataBosReads_passing_samples,
    output:
        'results/{library}/00_QC/seqkit.report.KDBos.txt'
    benchmark:
        'benchmarks/{library}/seqkitKneaddataHostReads.txt'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.HostReads} > {output} '


def get_seqkitKneaddataCapraReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/02_kneaddata/{samples}_CAPRA_ARS1.2_bowtie2_contam.fastq.gz", samples = passed, library = lib)


rule seqkitKneaddataCapraReads:
    input:
        HostReads = get_seqkitKneaddataCapraReads_passing_samples,
    output:
        'results/{library}/00_QC/seqkit.report.KDCapra.txt'
    benchmark:
        'benchmarks/{library}/seqkitKneaddataCapraReads.txt'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.HostReads} > {output} '

def get_seqkitKneaddataCervusReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/02_kneaddata/{samples}_mCerEla1_bowtie2_contam.fastq.gz", samples = passed, library = lib)


rule seqkitKneaddataCervusReads:
    input:
        HostReads = get_seqkitKneaddataCervusReads_passing_samples,
    output:
        'results/{library}/00_QC/seqkit.report.KDCervus.txt'
    benchmark:
        'benchmarks/{library}/seqkitKneaddataCervusReads.txt'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.HostReads} > {output} '


def get_seqkitKneaddataSILVAReads_passing_samples(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/02_kneaddata/{samples}_SLIVA138.1_bowtie2_contam.fastq.gz", samples = passed, library = lib)


rule seqkitKneaddataSILVAReads:
    input:
        silvaReads = get_seqkitKneaddataSILVAReads_passing_samples,
    output:
        'results/{library}/00_QC/seqkit.report.KDSILVA138.txt'
    benchmark:
        'benchmarks/{library}/seqkitKneaddataSILVAReads.txt'
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 30),
        partition='compute',
    shell:
        'seqkit stats -j {threads} -a {input.silvaReads} > {output} '


#KRAKEN2 RULES
rule kraken2GTDB:
    input:
        KDRs = "results/{library}/02_kneaddata/{samples}.fastq.gz",
    output:
        k2OutputGTDB = "results/{library}/03_kraken2GTDB/{samples}.k2",
        k2ReportGTDB = "results/{library}/03_kraken2GTDB/{samples}.kraken2",
    log:
        "logs/{library}/kraken2GTDB/kraken2GTDB.{samples}.GTDB.log",
    benchmark:
        "benchmarks/{library}/kraken2GTDB.{samples}.txt"
    conda:
        "kraken2"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 324 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
        partition = "compute,hugemem"
    shell:
        "kraken2 "
        "--use-names "
        "--db /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB "
        "-t {threads} "
        "--report {output.k2ReportGTDB} "
        "--report-minimizer-data "
        "{input.KDRs} "
        "--output {output.k2OutputGTDB} "
        " | tee {log} 2>&1 "


def get_prinseq_passing_samples_for_taxpasta(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.seqkitMaskingPrinseqReads.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand("results/{library}/03_kraken2GTDB/{samples}.kraken2", samples = passed, library = lib)

rule taxpastaKraken2:
    input:
        get_prinseq_passing_samples_for_taxpasta,
        #expand("results/{library}/03_kraken2GTDB/{samples}.kraken2", samples = FIDs, library = LIBRARY),
    output:
        "results/{library}/kraken2.{library}.genus.counts.tsv",
    benchmark:
        "benchmarks/{library}/taxpastaKraken2.txt"
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute,hugemem"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at genus "
        "{input} "
