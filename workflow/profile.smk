# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"

import os

FID, = glob_wildcards("results/02_kneaddata/{FID}.fastq")

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
        "results/centrifuge.counts.tsv",
        "results/centrifuge.counts.biom",

        "results/kraken2.counts.tsv",
        "results/kraken2.counts.biom",

        "results/bracken.k2.counts.tsv",
        "results/bracken.k2.counts.biom",

        expand("results/03_humann3Uniref50EC/{sample}_pathcoverage.tsv", sample=FID),

wildcard_constraints:
    sample="[a-zA-Z0-9]+"

localrules: generateCentrifugeSampleSheet


rule generateCentrifugeSampleSheet:
    output:
        sampleSheet='resources/centrifugeSampleSheet.tsv',
    threads:2
    shell: 
        './workflow/scripts/generate_centrifuge_sample_sheet.sh -d results/02_kneaddata -p kneaddata.trimmed.fastq -o {output.sampleSheet} '


rule centrifugeGTDB:
    input:
        sampleSheet = "resources/centrifugeSampleSheet.tsv",
    output:
        out = expand("results/03_centrifuge/{sample}.GTDB.centrifuge", sample = FID),
        report = expand("results/03_centrifuge/{sample}.GTDB.centrifuge.report", sample = FID),
    log:
        "logs/centrifuge.GTDB.multi.log",
    conda:
        "centrifuge"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 140 + ((attempt - 1) + 20),
        time = "06:00:00",
    shell:
        "centrifuge "
        "-x /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/centrifuge/GTDB " #TODO config
        "--sample-sheet {input.sampleSheet} "
        "-t "
        "--threads {threads} "
        "2>&1 | tee {log}"


rule centrifugeKrakenReport:
    input:
        centrifuge = "results/03_centrifuge/{sample}.GTDB.centrifuge",
    output:
        centrifugeKraken2 = "results/03_centrifuge/{sample}.centrifuge",
    log:
        "logs/centrifuge.to.kraken2.{sample}.log",
    conda:
        "centrifuge"
    threads: 2
    shell:
        "centrifuge-kreport "
        "-x /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/centrifuge/GTDB " #TODO config
        "{input.centrifuge} > "
        "{output.centrifugeKraken2}"


rule brackenCentrifugeGenus:
    input:
        centrifugeKraken2='results/03_centrifuge/{sample}.GTDB.centrifuge.k2report',
    output:
        braken='results/04_braken/{sample}.GTDB.centrifuge.k2report.T1.bracken.genus',
        brakenReport='results/04_braken/{sample}.GTDB.centrifuge.k2report.T1.bracken.genus.report',
    log:
        'logs/centrifuge.bracken.genus.{sample}.GTDB.log'
    conda:
        'kraken2'
    threads: 2 
    shell:
        'bracken '
        '-d /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/kraken/GTDB ' #TODO config
        '-i {input.centrifugeKraken2} '
        '-o {output.braken} '
        '-w {output.brakenReport} '
        '-r 80 '
        '-l G '
        '-t 1 '
        '&> {log} '


rule taxpastaCentrifugeTable:
    input:
        expand("results/03_centrifuge/{sample}.centrifuge", sample = FID),
    output:
        "results/centrifuge.counts.tsv",
    conda:
        "taxpasta"
    threads: 2
    shell:
        "taxpasta merge "
        "-p centrifuge "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/centrifuge " #TODO config
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at genus "
        "{input} "


rule taxpastaCentrifugeBiom:
    input:
        expand("results/03_centrifuge/{sample}.centrifuge", sample = FID),
    output:
        "results/centrifuge.counts.biom",
    conda:
        "taxpasta"
    threads: 2
    shell:
        "taxpasta merge "
        "-p centrifuge "
        "-o {output} "
        "--output-format BIOM "
        "--taxonomy /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/centrifuge " #TODO config
        "--add-name "
        "--summarise-at genus "
        "{input} "


rule brackenCentrifugeSpecies:
    input:
        centrifugeKraken2='results/03_centrifuge/{sample}.GTDB.centrifuge.k2report',
    output:
        braken='results/04_braken/{sample}.GTDB.centrifuge.k2report.T1.bracken.species',
        brakenReport='results/04_braken/{sample}.GTDB.centrifuge.k2report.T1.bracken.species.report',
    log:
        'logs/centrifuge.bracken.species.{samples}.GTDB.log'
    conda:
        'kraken2'
    threads: 2 
    shell:
        'bracken '
        '-d /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/kraken/GTDB ' #TODO config
        '-i {input.centrifugeKraken2} '
        '-o {output.braken} '
        '-w {output.brakenReport} '
        '-r 80 '
        '-l S '
        '-t 1 '
        '&> {log} '


#KRAKEN2 RULES
rule kraken2GTDB:
    input:
        KDRs = "results/02_kneaddata/{sample}.fastq",
    output:
        k2OutputGTDB = "results/03_kraken2GTDB/{sample}.k2",
        k2ReportGTDB = "results/03_kraken2GTDB/{sample}.kraken2",
    log:
        "logs/kraken2GTDB.{sample}.GTDB.log",
    conda:
        "kraken2"
    threads: 18
    resources:
        # dynamic memory allocation: start with 400G and increment by 20G with every failed attempt 
        mem_gb = lambda wildcards, attempt: 360 + ((attempt - 1) * 20),
    shell:
        "kraken2 "
        "--use-names "
        "--quick "
        "--db /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB " #TODO config
        "-t {threads} "
        "--report {output.k2ReportGTDB} "
        "--report-minimizer-data "
        "{input.KDRs} "
        "--output {output.k2OutputGTDB} "
        "2>&1 | tee {log} "


rule taxpastaKraken2:
    input:
        expand("results/03_kraken2GTDB/{sample}.kraken2", sample = FID),
    output:
        "results/kraken2.counts.tsv",
    conda:
        "taxpasta"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB/taxonomy " #TODO config
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at genus "
        "{input} "


rule taxpastaKraken2Biom:
    input:
        expand("results/03_kraken2GTDB/{sample}.kraken2", sample = FID),
    output:
        "results/kraken2.counts.biom",
    conda:
        "taxpasta"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format BIOM "
        "--taxonomy /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB/taxonomy " #TODO config
        "--add-name "
        "--summarise-at genus "
        "{input} "

# BRACKEN RULES
rule brackenSpecies:
    input:
        k2ReportGTDB = "results/03_kraken2GTDB/{sample}.kraken2",
    output:
        bOutput = "results/03_brackenSpecies/{sample}.bracken",
        bReport = "results/03_brackenSpecies/{sample}.br",
    log:
        "logs/brackenSpecies.{sample}.log",
    conda:
        "kraken2"
    threads: 2
    shell:
        "bracken "
        "-d /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB " #TODO config
        "-i {input.k2ReportGTDB} "
        "-o {output.bOutput} "
        "-w {output.bReport} "
        "-r 80 "
        "-l S "
        "-t 10 " # Necessary to get floating point counts to sum to 1.0 in taxpasta
        "&> {log} "


rule taxpastaKraken2Bracken:
    input:
        expand("results/03_brackenSpecies/{sample}.bracken", sample = FID),
    output:
        "results/bracken.k2.counts.tsv",
    conda:
        "taxpasta"
    shell:
        "taxpasta merge "
        "-p bracken "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB/taxonomy " #TODO config
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at species "
        "{input} "


rule taxpastaKraken2BrackenBiom:
    input:
        expand("results/03_brackenSpecies/{sample}.bracken", sample = FID),
    output:
        "results/bracken.k2.counts.biom",
    conda:
        "taxpasta"
    shell:
        "taxpasta merge "
        "-p bracken "
        "-o {output} "
        "--output-format BIOM "
        "--taxonomy /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB/taxonomy " #TODO config
        "--add-name "
        "--summarise-at species "
        "{input} "


# HUMANN# RULES
rule humann3Uniref50EC:
    input:
        kneaddataReads = "results/02_kneaddata/{sample}.fastq",
    output:
        genes = "results/03_humann3Uniref50EC/{sample}_genefamilies.tsv",
        pathways = "results/03_humann3Uniref50EC/{sample}_pathabundance.tsv",
        pathwaysCoverage = "results/03_humann3Uniref50EC/{sample}_pathcoverage.tsv",
    log:
        "logs/humann3.{sample}.uniref50EC.log",
    conda:
        "biobakery"
    threads: 16
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) + 12),
    message:
        "humann3 profiling with uniref50EC: {wildcards.samples}\n"
    shell:
        "humann3 "
        "--memory-use maximum "
        "--threads {threads} "
        "--bypass-nucleotide-search "
        "--search-mode uniref50 "
        "--protein-database /bifo/scratch/2022-BJP-GTDB/biobakery/humann3/unirefECFilt " #TODO config
        "--input-format fastq "
        "--output results/03_humann3Uniref50EC "
        "--input {input.kneaddataReads} "
        "--output-basename {wildcards.samples} "
        "--o-log {log} "
        "--remove-temp-output "
