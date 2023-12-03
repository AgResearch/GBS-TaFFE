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
kdr_seqkit_report = os.path.join("results", LIBRARY, "00_QC", "seqkit.report.KDR.txt")


def get_passing_KDR_files_GTDB207(wildcards, seqkitOut = kdr_seqkit_report, minReads=min_reads, lib=LIBRARY):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "03_kraken2GTDB/{samples}.GTDB207.kraken2"), samples = passed)

def get_passing_KDR_files_GTDB214(wildcards, seqkitOut = kdr_seqkit_report, minReads=min_reads, lib=LIBRARY):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "03_kraken2GTDB/{samples}.GTDB214.kraken2"), samples = passed)

def get_passing_FIDs(seqkitOut = kdr_seqkit_report, minReads=min_reads, lib=LIBRARY):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    return qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()

FIDs = get_passing_FIDs()

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
        #GTDB207
        # expand("results/{library}/kraken2.GTDB207.kingdom.counts.tsv",  library = LIBRARY),
        # expand("results/{library}/kraken2.GTDB207.phylum.counts.tsv",  library = LIBRARY),
        # expand("results/{library}/kraken2.GTDB207.order.counts.tsv",  library = LIBRARY),
        # expand("results/{library}/kraken2.GTDB207.class.counts.tsv",  library = LIBRARY),
        # expand("results/{library}/kraken2.GTDB207.family.counts.tsv",  library = LIBRARY),
        # expand("results/{library}/kraken2.GTDB207.genus.counts.tsv",  library = LIBRARY),
        # expand("results/{library}/kraken2.GTDB207.species.counts.tsv",  library = LIBRARY),
        # expand("results/{library}/kraken2.GTDB207.genus.counts.biom",  library = LIBRARY),
        #GTDB214
        expand("results/{library}/kraken2.GTDB214.kingdom.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.phylum.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.order.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.class.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.family.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.genus.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.species.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.genus.counts.biom",  library = LIBRARY),


        # "results/04_functional/humann3_uniref50EC_pathabundance.rpk.tsv",
        # "results/04_functional/humann3_uniref50EC_genefamilies.rpk.tsv",
        # "results/04_functional/humann3_uniref50EC_pathcoverage.rpk.tsv",
        # "results/04_functional/humann3_uniref50EC_genefamilies.rpk.KO.tsv",
        # "results/04_functional/humann3_uniref50EC_genefamilies.rpk.EC.tsv",
        # "results/04_functional/humann3_uniref50EC_genefamilies.rpk.pfam.tsv",
        # "results/04_functional/humann3_uniref50EC_genefamilies.rpk.EggNOG.tsv",
        # "results/04_functional/humann3_uniref50EC_pathabundance.rpk.cpm.QC.tsv",
        # "results/04_functional/humann3_uniref50EC_genefamilies.rpk.cpm.QC.tsv"

        #expand("results/03_humann3Uniref50EC/{sample}_pathcoverage.tsv", sample=FIDs),


#KRAKEN2 RULES
rule kraken2_GTDB207:
    input:
        KDRs = "results/{library}/02_kneaddata/{samples}.fastq.gz",
    output:
        k2OutputGTDB = "results/{library}/03_kraken2GTDB/{samples}.GTDB207.k2",
        k2ReportGTDB = "results/{library}/03_kraken2GTDB/{samples}.GTDB207.kraken2",
    log:
        "logs/{library}/kraken2_GTDB207/kraken2GTDB.{samples}.GTDB207.log",
    benchmark:
        "benchmarks/{library}/kraken2_GTDB207.{samples}.txt"
    conda:
        "kraken2"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 324 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
        partition = "hugemem"
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



rule taxpasta_Kraken2_GTDB207_kingdom:
    input:
        get_passing_KDR_files_GTDB207,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB207.kingdom.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB207_kingdom.txt")
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
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at kingdom "
        "{input} "


rule taxpasta_Kraken2_GTDB207_phylum:
    input:
        get_passing_KDR_files_GTDB207,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB207.phylum.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB207_phylum.txt")
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
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at phylum "
        "{input} "




rule taxpasta_Kraken2_GTDB207_order:
    input:
        get_passing_KDR_files_GTDB207,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB207.order.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB207_order.txt")
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
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at order "
        "{input} "


rule taxpasta_Kraken2_GTDB207_class:
    input:
        get_passing_KDR_files_GTDB207,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB207.class.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB207_class.txt")
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
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at class "
        "{input} "


rule taxpasta_Kraken2_GTDB207_family:
    input:
        get_passing_KDR_files_GTDB207,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB207.family.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB207_family.txt")
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
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at family "
        "{input} "


rule taxpasta_Kraken2_GTDB207_genus:
    input:
        get_passing_KDR_files_GTDB207,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB207.genus.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB207_genus.txt")
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
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at genus "
        "{input} "


rule taxpasta_Kraken2_GTDB207_species:
    input:
        get_passing_KDR_files_GTDB207,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB207.species.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB207_species.txt")
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
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at species "
        "{input} "


rule taxpasta_Kraken2_GTDB207_Biom:
    input:
        get_passing_KDR_files_GTDB207,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB207.genus.counts.biom")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB207_Biom.txt")
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
        "--taxonomy /agr/scratch/projects/2023-mbie-rumen-gbs/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--summarise-at genus "
        "{input} "


#KRAKEN2 RULES
rule kraken2_GTDB214:
    input:
        KDRs = "results/{library}/02_kneaddata/{samples}.fastq.gz",
    output:
        k2OutputGTDB = "results/{library}/03_kraken2GTDB/{samples}.GTDB214.k2",
        k2ReportGTDB = "results/{library}/03_kraken2GTDB/{samples}.GTDB214.kraken2",
    log:
        "logs/{library}/kraken2_GTDB214/kraken2GTDB.{samples}.GTDB214.log",
    benchmark:
        "benchmarks/{library}/kraken2_GTDB214.{samples}.txt"
    conda:
        "kraken2"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 324 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 30),
        partition = "hugemem"
    shell:
        "kraken2 "
        "--use-names "
        "--db /agr/scratch/projects/2022-bjp-gtdb/build-GTDB-DBs/GTDB/kraken2-GTDB-214.1 " 
        "-t {threads} "
        "--report {output.k2ReportGTDB} "
        "--report-minimizer-data "
        "{input.KDRs} "
        "--output {output.k2OutputGTDB} "
        " | tee {log} 2>&1 "



rule taxpasta_Kraken2_GTDB214_kingdom:
    input:
        get_passing_KDR_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.kingdom.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB214_kingdom.txt")
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
        "--summarise-at kingdom "
        "{input} "


rule taxpasta_Kraken2_GTDB214_phylum:
    input:
        get_passing_KDR_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.phylum.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB214_phylum.txt")
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
    input:
        get_passing_KDR_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.order.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB214_order.txt")
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
    input:
        get_passing_KDR_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.class.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB214_class.txt")
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
    input:
        get_passing_KDR_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.family.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB214_family.txt")
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
    input:
        get_passing_KDR_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.genus.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB214_genus.txt")
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
    input:
        get_passing_KDR_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.species.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB214_species.txt")
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
    input:
        get_passing_KDR_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.genus.counts.biom")
    benchmark:
        os.path.join("results", LIBRARY, "taxpasta_Kraken2_GTDB214_Biom.txt")
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



