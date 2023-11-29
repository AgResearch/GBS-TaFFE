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


def get_passing_KDR_files(wildcards, seqkitOut = "results/00_QC/seqkit.report.KDR.txt", minReads=min_reads, lib=LIBRARY):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "03_kraken2GTDB/{samples}.kraken2"), samples = passed)


def get_passing_FIDs(wildcards, seqkitOut = "results/00_QC/seqkit.report.KDR.txt", minReads=min_reads, lib=LIBRARY):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    return qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()

FIDs = get_passing_FIDs(seqkitOut = "results/00_QC/seqkit.report.KDR.txt", minReads = min_reads, lib = LIBRARY)

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
        expand("results/{library}/kraken2.genus.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.genus.counts.biom",  library = LIBRARY),
        #"results/bracken.k2.counts.tsv",
        #"results/bracken.k2.counts.biom",
        #expand("results/03_centrifuge/{sample}.centrifuge", sample=FIDs),
        #expand("results/03_humann3Uniref50EC/{sample}_pathcoverage.tsv", sample=FIDs),


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
        "--db /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB " 
        "-t {threads} "
        "--report {output.k2ReportGTDB} "
        "--report-minimizer-data "
        "{input.KDRs} "
        "--output {output.k2OutputGTDB} "
        " | tee {log} 2>&1 "



rule taxpastaKraken2:
    input:
        get_passing_KDR_files,
    output:
        os.path.join("results", LIBRARY, "kraken2.genus.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpastaKraken2.txt")
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
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at genus "
        "{input} "


rule taxpastaKraken2Biom:
    input:
        get_passing_KDR_files,
    output:
        os.path.join("results", LIBRARY, "kraken2.genus.counts.biom")
    benchmark:
        os.path.join("results", LIBRARY, "taxpastaKraken2Biom.txt")
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
        "--output-format BIOM "
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--summarise-at genus "
        "{input} "


# # HUMANN RULES #TODO Import code from TaFFE for aggregating and converting the outputs
# rule humann3Uniref50EC:
#     input:
#         kneaddataReads = "results/02_kneaddata/{sample}.fastq.gz",
#     output:
#         genes = "results/03_humann3Uniref50EC/{sample}_genefamilies.tsv",
#         pathways = "results/03_humann3Uniref50EC/{sample}_pathabundance.tsv",
#         pathwaysCoverage = "results/03_humann3Uniref50EC/{sample}_pathcoverage.tsv",
#     log:
#         "logs/humann3.{sample}.uniref50EC.log",
#     benchmark:
#         "benchmarks/humann3Uniref50EC.{sample}.txt"
#     conda:
#         "biobakery"
#     threads: 16
#     resources:
#         mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) + 12),
#         time = lambda wildcards, attempt: 24 + ((attempt - 1) + 12),
#         partition = "large,milan"
#     message:
#         "humann3 profiling with uniref50EC: {wildcards.sample}\n"
#     shell:
#         "humann3 "
#         "--memory-use maximum "
#         "--threads {threads} "
#         "--bypass-nucleotide-search "
#         "--search-mode uniref50 "
#         "--protein-database /nesi/nobackup/agresearch03843/biobakery/biobakery/humann3/unirefECFilt "
#         "--input-format fastq "
#         "--output results/03_humann3Uniref50EC "
#         "--input {input.kneaddataReads} "
#         "--output-basename {wildcards.sample} "
#         "--o-log {log} "
#         "--remove-temp-output "



# BRACKEN RULES #TODO Implement
# rule brackenSpecies:
#     input:
#         k2ReportGTDB = "results/03_kraken2GTDB/{sample}.kraken2",
#     output:
#         bOutput = "results/03_brackenSpecies/{sample}.bracken",
#         bReport = "results/03_brackenSpecies/{sample}.br",
#     log:
#         "logs/brackenSpecies.{sample}.log",
#     benchmark:
#         "benchmarks/brackenSpecies.{sample}.txt"
#     conda:
#         "kraken2"
#     threads: 2
#     resources:
#         mem_gb = 1,
#         time = 2,
#         partition = "large,milan"
#     shell:
#         "bracken "
#         "-d /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB "
#         "-i {input.k2ReportGTDB} "
#         "-o {output.bOutput} "
#         "-w {output.bReport} "
#         "-r 80 "
#         "-l S "
#         # "-t 10 "
#         "&> {log} "


# rule taxpastaKraken2Bracken:
#     input:
#         expand("results/03_brackenSpecies/{sample}.bracken", sample = FIDs),
#     output:
#         "results/bracken.k2.counts.tsv",
#     benchmark:
#         "benchmarks/taxpastaKraken2Bracken.txt"
#     conda:
#         "taxpasta"
#     threads: 2
#     resources:
#         mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
#         time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
#         partition = "large,milan"
#     shell:
#         "taxpasta merge "
#         "-p bracken "
#         "-o {output} "
#         "--output-format TSV "
#         "--taxonomy /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB/taxonomy "
#         "--add-name "
#         "--add-rank "
#         "--add-lineage "
#         "--summarise-at species "
#         "{input} "


# rule taxpastaKraken2BrackenBiom:
#     input:
#         expand("results/03_brackenSpecies/{sample}.bracken", sample = FIDs),
#     output:
#         "results/bracken.k2.counts.biom",
#     benchmark:
#         "benchmarks/taxpastaKraken2BrackenBiom.txt"
#     conda:
#         "taxpasta"
#     threads: 2
#     resources:
#         mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
#         time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
#         partition = "large,milan"
#     shell:
#         "taxpasta merge "
#         "-p bracken "
#         "-o {output} "
#         "--output-format BIOM "
#         "--taxonomy /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB/taxonomy "
#         "--add-name "
#         "--summarise-at species "
#         "{input} "



