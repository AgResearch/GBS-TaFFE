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
seqkit_report = os.path.join("results", LIBRARY, "00_QC", "seqkit.report.raw.txt")

def get_passing_files_GTDB214(wildcards, seqkitOut = seqkit_report, minReads=min_reads, lib=LIBRARY):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > minReads]
    passed = qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()
    return expand(os.path.join("results", lib, "03_kraken2_GTDB214/{samples}.GTDB214.kraken2"), samples = passed)

def get_passing_FIDs(seqkitOut = seqkit_report, minReads=min_reads, lib=LIBRARY):
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
#kraken2
        expand("results/{library}/kraken2.GTDB214.domain.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.phylum.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.order.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.class.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.family.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.genus.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.species.counts.tsv",  library = LIBRARY),
        expand("results/{library}/kraken2.GTDB214.genus.counts.biom",  library = LIBRARY),
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


#KRAKEN2 RULES
rule kraken2_GTDB214:
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
    priority: 1000
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


rule taxpasta_Kraken2_GTDB214_domain:
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.domain.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "taxpasta_Kraken2_GTDB214_kingdom.txt")
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
    input:
        get_passing_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.phylum.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "taxpasta_Kraken2_GTDB214_phylum.txt")
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
        get_passing_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.order.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "taxpasta_Kraken2_GTDB214_order.txt")
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
        get_passing_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.class.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "taxpasta_Kraken2_GTDB214_class.txt")
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
        get_passing_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.family.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "taxpasta_Kraken2_GTDB214_family.txt")
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
        get_passing_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.genus.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "taxpasta_Kraken2_GTDB214_genus.txt")
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
        get_passing_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.species.counts.tsv")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "taxpasta_Kraken2_GTDB214_species.txt")
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
        get_passing_files_GTDB214,
    output:
        os.path.join("results", LIBRARY, "kraken2.GTDB214.genus.counts.biom")
    benchmark:
        os.path.join("results", LIBRARY, "benchmarks", "taxpasta_Kraken2_GTDB214_Biom.txt")
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
    priority: 1000
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


def get_humann3_pathcoverage(wildcards, seqkitOut = seqkit_report, minReads=min_reads, lib=LIBRARY):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
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


def get_humann3_genefamilies(wildcards, seqkitOut = seqkit_report, minReads=min_reads, lib=LIBRARY):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
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


def get_humann3_pathcoverage(wildcards, seqkitOut = seqkit_report, minReads=min_reads, lib=LIBRARY):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
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


