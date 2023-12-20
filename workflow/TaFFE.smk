# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"

import os
import pandas as pd

wildcard_constraints:
    samples="\w+"


# Global minimum read count for processing
min_reads = 25000

LIBRARY = config["LIBRARY"]

# seqkit_report = os.path.join("results", LIBRARY, "00_QC", "seqkit.report.raw.txt")

# def get_passing_FIDs(seqkitOut = seqkit_report, minReads=min_reads, lib=LIBRARY):
#     import pandas as pd
#     qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
#     qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
#     qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
#     return qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()

# FIDs = get_passing_FIDs()


input_fastq_pattern = os.path.join('results', config["LIBRARY"], '01_cutadapt', '{samples}.fastq.gz')
print(input_fastq_pattern)
FIDs, = glob_wildcards(input_fastq_pattern)


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


rule all:
    input:
        expand('results/{library}/00_QC/seqkit.report.KDR.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.raw.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.bbduk.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.prinseq.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDTrim.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDTRF.txt', library = LIBRARY),
        expand('results/{library}/00_QC/seqkit.report.KDSILVA138.txt', library = LIBRARY),
#kraken2
        expand("results/{library}/{library}.kraken2.GTDB214.domain.counts.tsv",  library = LIBRARY),
        expand("results/{library}/{library}.kraken2.GTDB214.phylum.counts.tsv",  library = LIBRARY),
        expand("results/{library}/{library}.kraken2.GTDB214.order.counts.tsv",  library = LIBRARY),
        expand("results/{library}/{library}.kraken2.GTDB214.class.counts.tsv",  library = LIBRARY),
        expand("results/{library}/{library}.kraken2.GTDB214.family.counts.tsv",  library = LIBRARY),
        expand("results/{library}/{library}.kraken2.GTDB214.genus.counts.tsv",  library = LIBRARY),
        # expand("results/{library}/{library}.kraken2.GTDB214.species.counts.tsv",  library = LIBRARY),
        expand("results/{library}/{library}.kraken2.GTDB214.genus.counts.biom",  library = LIBRARY),
#human3
        # expand("results/{library}/05_functional/humann3_uniref50EC_microbial_pathabundance.rpk.tsv", library = LIBRARY),
        # expand("results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.tsv", library = LIBRARY),
        # expand("results/{library}/05_functional/humann3_uniref50EC_microbial_pathcoverage.rpk.tsv", library = LIBRARY),
        # expand("results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.KO.tsv", library = LIBRARY),
        # expand("results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.EC.tsv", library = LIBRARY),
        # expand("results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.pfam.tsv", library = LIBRARY),
        # expand("results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.EggNOG.tsv", library = LIBRARY),
        # expand("results/{library}/05_functional/humann3_uniref50EC_microbial_pathabundance.rpk.cpm.QC.tsv", library = LIBRARY),
        # expand("results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.cpm.QC.tsv", library = LIBRARY),


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
    priority: 10000
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
    priority: 10000
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
    priority: 10000
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
    priority: 10000
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
    priority: 10000
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
    priority: 10000
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


#KRAKEN2 RULES
rule kraken2_GTDB214:
    priority: 1000
    input:
        KDRs = "results/{library}/02_kneaddata/{samples}.fastq.gz",
    output:
        k2OutputGTDB = "results/{library}/03_kraken2_GTDB214/{samples}.GTDB214.k2",
        k2ReportGTDB = "results/{library}/03_kraken2_GTDB214/{samples}.GTDB214.kraken2",
        k2Classified_read = temp("results/{library}/03_kraken2_GTDB214/{samples}.GTDB214.kraken2.classified.fastq"),
    log:
        os.path.join("results", "{library}", "logs", "kraken2_GTDB214", "kraken2GTDB.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results","{library}", "benchmarks", "kraken2_GTDB214.{samples}.txt"),
    conda:
        "kraken2"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 420 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
        partition = "hugemem,compute"
    shell:
        "kraken2 "
        "--use-names "
        "--report-zero-counts "
        "--gzip-compressed "
        "--classified-out {output.k2Classified_read} "
        "--db /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1 " 
        "-t {threads} "
        "--report {output.k2ReportGTDB} "
        "--report-minimizer-data "
        "--output {output.k2OutputGTDB} "
        "{input.KDRs} "
        "2>&1 | tee {log} "


rule kraken2_GTDB214_gz:
    input:
        k2OutputGTDB = "results/{library}/03_kraken2_GTDB214/{samples}.GTDB214.k2",
        k2Classified_read = "results/{library}/03_kraken2_GTDB214/{samples}.GTDB214.kraken2.classified.fastq",
    output:
        k2OutputGTDB = "results/{library}/03_kraken2_GTDB214/{samples}.GTDB214.k2.gz",
        k2Classified_read = temp("results/{library}/03_kraken2_GTDB214/{samples}.GTDB214.kraken2.classified.fastq.gz"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_GTDB214_gz.{samples}.txt"),
    conda:
        "pigz"
    threads: 16
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 16),
        time = lambda wildcards, attempt: 40 + ((attempt - 1) * 60),
        partition = "compute,hugemem"
    shell:
        "pigz "
        "-p {threads} "
        "{input.k2OutputGTDB} " 
        "{input.k2Classified_read} " 


def get_passing_files_GTDB214(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "03_kraken2_GTDB214/{samples}.GTDB214.kraken2"), samples = passed)


# def get_passing_files_GTDB214(wildcards, seqkitOut = seqkit_report, minReads=min_reads, lib=LIBRARY):
#     import pandas as pd
#     qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
#     qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
#     qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
#     passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
#     return expand(os.path.join("results", lib, "03_kraken2_GTDB214/{samples}.GTDB214.kraken2"), samples = passed)


rule taxpasta_Kraken2_GTDB214_domain:
    priority: 1000
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", "{library}", "{library}.kraken2.GTDB214.domain.counts.tsv")
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "taxpasta_Kraken2_GTDB214_kingdom.txt")
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at domain "
        "{input} "


rule taxpasta_Kraken2_GTDB214_phylum:
    priority: 1000
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", "{library}", "{library}.kraken2.GTDB214.phylum.counts.tsv")
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "taxpasta_Kraken2_GTDB214_phylum.txt")
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at phylum "
        "{input} "


rule taxpasta_Kraken2_GTDB214_order:
    priority: 1000
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", "{library}", "{library}.kraken2.GTDB214.order.counts.tsv")
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "taxpasta_Kraken2_GTDB214_order.txt")
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at order "
        "{input} "


rule taxpasta_Kraken2_GTDB214_class:
    priority: 1000
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", "{library}", "{library}.kraken2.GTDB214.class.counts.tsv")
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "taxpasta_Kraken2_GTDB214_class.txt")
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at class "
        "{input} "


rule taxpasta_Kraken2_GTDB214_family:
    priority: 1000
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", "{library}", "{library}.kraken2.GTDB214.family.counts.tsv")
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "taxpasta_Kraken2_GTDB214_family.txt")
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at family "
        "{input} "


rule taxpasta_Kraken2_GTDB214_genus:
    priority: 1000
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", "{library}", "{library}.kraken2.GTDB214.genus.counts.tsv")
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "taxpasta_Kraken2_GTDB214_genus.txt")
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at genus "
        "{input} "


rule taxpasta_Kraken2_GTDB214_species:
    priority: 1000
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", "{library}", "{library}.kraken2.GTDB214.species.counts.tsv")
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "taxpasta_Kraken2_GTDB214_species.txt")
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at species "
        "{input} "


rule taxpasta_Kraken2_GTDB214_Biom:
    priority: 1000
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", "{library}", "{library}.kraken2.GTDB214.genus.counts.biom")
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "taxpasta_Kraken2_GTDB214_Biom.txt")
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format BIOM "
        "--taxonomy /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1/taxonomy "
        "--add-name "
        "--summarise-at genus "
        "{input} "


# HUMANN RULES
rule kraken2_host_filter:
    input:
        k2Classified_read = "results/{library}/03_kraken2_GTDB214/{samples}.GTDB214.kraken2.classified.fastq.gz",
    output:
        k2OutputHosts = temp("results/{library}/04_k2_filtering/{samples}.Hosts.k2"),
        k2_filtered_read = temp("results/{library}/04_k2_filtering/{samples}.nonhost.fastq"),
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter.{samples}.txt"),
    conda:
        "kraken2"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 32),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
        partition = "compute,hugemem"
    shell:
        "kraken2 "
        "--gzip-compressed "
        "--unclassified-out {output.k2_filtered_read} "
        "--db /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-hosts " 
        "-t {threads} "
        "--output {output.k2OutputHosts} "
        "{input.k2Classified_read} "
        "2>&1 | tee {log} "


rule kraken2_host_filter_gz:
    input:
        k2_filtered_read = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq",
    output:
        k2_filtered_reads_gz = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq.gz",
    log:
        os.path.join("results", "{library}", "logs", "kraken2", "kraken2_host_filter_gz.{samples}.GTDB214.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "kraken2_host_filter_gz.{samples}.txt"),
    conda:
        "pigz"
    threads: 16
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 16),
        time = lambda wildcards, attempt: 20 + ((attempt - 1) * 20),
        partition = "compute,hugemem"
    shell:
        "pigz -p 12 "
        "{input.k2_filtered_read} "


rule humann3Uniref50EC:
    input:
        k2_filtered_reads_gz = "results/{library}/04_k2_filtering/{samples}.nonhost.fastq.gz",
    output:
        genes = "results/{library}/05_humann3Uniref50EC/{samples}.genefamilies.tsv",
        pathways = "results/{library}/05_humann3Uniref50EC/{samples}.pathabundance.tsv",
        pathwaysCoverage = "results/{library}/05_humann3Uniref50EC/{samples}.pathcoverage.tsv",
    log:
        os.path.join("results", "{library}", "logs", "humann3", "humann3.Uniref50EC.{samples}.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "humann3.Uniref50EC.{samples}.txt"),
    conda:
        "humann3"
    threads: 24
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 480 + ((attempt - 1) + 480),
        partition = "compute,hugemem"
    shell:
        "humann3 "
        "--memory-use maximum "
        "--threads {threads} "
        "--bypass-nucleotide-search "
        "--search-mode uniref50 "
        "--protein-database /agr/scratch/projects/2023-mbie-rumen-gbs/biobakery/biobakery/humann3/uniref50ECFilt "
        "--input-format fastq.gz "
        "--output results/{LIBRARY}/05_humann3Uniref50EC "
        "--input {input.k2_filtered_reads_gz} "
        "--output-basename {wildcards.samples} "
        "--o-log {log} "
        "--remove-temp-output && "
        "mv results/{LIBRARY}/05_humann3Uniref50EC/{wildcards.samples}_genefamilies.tsv {output.genes}; "
        "mv results/{LIBRARY}/05_humann3Uniref50EC/{wildcards.samples}_pathabundance.tsv {output.pathways}; "
        "mv results/{LIBRARY}/05_humann3Uniref50EC/{wildcards.samples}_pathcoverage.tsv {output.pathwaysCoverage}; "


def get_humann3_pathcoverage(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "05_humann3Uniref50EC", "{samples}.pathcoverage.tsv"), samples = passed)

rule merge_functional_profiles_pathabundance:
    input:
        get_humann3_pathcoverage,
    output:
        os.path.join("results", LIBRARY, "05_functional", "humann3_uniref50EC_microbial_pathabundance.rpk.tsv")
    log:
        os.path.join("results", LIBRARY, "logs", "humann3.uniref50EC.merge.pathabundance.log")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "humann3.uniref50EC.merge.pathabundance.txt")
    conda:
        "humann3"
    threads:2
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) + 60),
        partition = "compute,hugemem"
    shell:
        "humann_join_tables "
        "-v "
        "-i results/{LIBRARY}/05_humann3Uniref50EC "
        "--file_name pathabundance "
        "-o {output} "
        "2>&1 | tee {log}  "


def get_humann3_genefamilies(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "05_humann3Uniref50EC", "{samples}.genefamilies.tsv"), samples = passed)

rule merge_functional_profiles_genefamilies:
    input:
        get_humann3_genefamilies
    output:
        os.path.join("results", LIBRARY, "05_functional", "humann3_uniref50EC_microbial_genefamilies.rpk.tsv")
    log:
        os.path.join("results", LIBRARY, "logs", "humann3.uniref50EC.merge.genefamilies.log")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "humann3.uniref50EC.merge.genefamilies.txt")
    conda:
        "humann3"
    threads:2
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) + 60),
        partition = "compute,hugemem"
    shell:
        "humann_join_tables "
        "-v "
        "-i results/{LIBRARY}/05_humann3Uniref50EC "
        "--file_name genefamilies "
        "-o {output} "
        "2>&1 | tee {log}  "


def get_humann3_pathcoverage(wildcards, minReads=min_reads, lib=LIBRARY):
    file = checkpoints.report_seqkit_raw.get().output[0]
    qc_stats = pd.read_csv(file, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "05_humann3Uniref50EC", "{samples}.pathcoverage.tsv"), samples = passed)

rule merge_functional_profiles_pathcoverage:
    input:
        get_humann3_pathcoverage,
    output:
        os.path.join("results", LIBRARY, "05_functional", "humann3_uniref50EC_microbial_pathcoverage.rpk.tsv")
    log:
        os.path.join("results", LIBRARY, "logs", "humann3", "humann3.uniref50EC.merge.pathcoverage.log")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "humann3.uniref50EC.merge.pathcoverage.txt")
    conda:
        "humann3"
    threads:2
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) + 60),
        partition = "compute,hugemem"
    shell:
        "humann_join_tables "
        "-v "
        "-i results/{LIBRARY}/05_humann3Uniref50EC "
        "--file_name pathcoverage "
        "-o {output} "
        "2>&1 | tee {log}  "


rule regroup_table_KO:
    input:
        genefamilies = "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.tsv"
    output:
        "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.KO.tsv"
    log:
        os.path.join("results", "{library}", "logs", "humann3", "humann3.uniref50EC.regroup.KO.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "humann3.uniref50EC.regroup.KO.log"),
    conda:
        "humann3"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) + 60),
        partition = "compute,hugemem"
    shell:
        "humann_regroup_table "
        "-i {input.genefamilies} "
        "--groups uniref50_ko "
        "--function sum "
        "--ungrouped Y "
        "--protected Y "
        "-o {output} "


rule regroup_table_EC:
    input:
        genefamilies = "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.tsv"
    output:
        "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.EC.tsv"
    log:
        os.path.join("results", "{library}", "logs", "humann3", "humann3.uniref50EC.regroup.EC.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "humann3.uniref50EC.regroup.EC.log"),
    conda:
        "humann3"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) + 60),
        partition = "compute,hugemem"
    shell:
        "humann_regroup_table "
        "-i {input.genefamilies} "
        "--groups uniref50_level4ec "
        "--function sum "
        "--ungrouped Y "
        "--protected Y "
        "-o {output} "


rule regroup_table_pfam:
    input:
        genefamilies = "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.tsv"
    output:
        "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.pfam.tsv"
    log:
        os.path.join("results", "{library}", "logs", "humann3", "humann3.uniref50EC.regroup.pfam.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "humann3.uniref50EC.regroup.pfam.log"),
    conda:
        "humann3"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) + 60),
        partition = "compute,hugemem"
    shell:
        "humann_regroup_table "
        "-i {input.genefamilies} "
        "--groups uniref50_pfam "
        "--function sum "
        "--ungrouped Y "
        "--protected Y "
        "-o {output} "


rule regroup_table_EggNOG:
    input:
        genefamilies = "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.tsv"
    output:
        "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.EggNOG.tsv"
    log:
        os.path.join("results", "{library}", "logs", "humann3", "humann3.uniref50EC.regroup.EggNOG.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "humann3.uniref50EC.regroup.EggNOG.txt"),
    conda:
        "humann3"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) + 60),
        partition = "compute,hugemem"
    shell:
        "humann_regroup_table "
        "-i {input.genefamilies} "
        "--groups uniref50_eggnog "
        "--function sum "
        "--ungrouped Y "
        "--protected Y "
        "-o {output} "


rule norm_humann3_genefamilies:
    input:
        genefamilies = "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.tsv"
    output:
        "results/{library}/05_functional/humann3_uniref50EC_microbial_genefamilies.rpk.cpm.QC.tsv"
    log:
        os.path.join("results", "{library}", "logs", "humann3", "humann3.uniref50EC.norm.genefamilies.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "humann3.uniref50EC.norm.genefamilies.txt"),
    conda:
        "humann3"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) + 60),
        partition = "compute,hugemem"
    shell:
        "humann_renorm_table "
        "-i {input.genefamilies} "
        "-u cpm "
        "-p "
        "-o {output} "


rule norm_humann3_pathabundance:
    input:
        pathabundance = "results/{library}/05_functional/humann3_uniref50EC_microbial_pathabundance.rpk.tsv"
    output:
        "results/{library}/05_functional/humann3_uniref50EC_microbial_pathabundance.rpk.cpm.QC.tsv"
    log:
        os.path.join("results", "{library}", "logs", "humann3", "humann3.uniref50EC.norm.pathabundance.log"),
    benchmark:
        os.path.join("results", "{library}", "benchmarks", "humann3.uniref50EC.norm.pathabundance.txt"),
    conda:
        "humann3"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) + 24),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) + 60),
        partition = "compute,hugemem"
    shell:
        "humann_renorm_table "
        "-i {input.pathabundance} "
        "-u cpm "
        "-p "
        "-o {output} "


