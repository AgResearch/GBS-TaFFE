# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"

import os

FID, = glob_wildcards("results/02_kneaddata/{fid}_kneaddata.trimmed.fastq")

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
        'results/00_qc/seqkit.report.raw.txt',
        'results/00_qc/seqkit.report.filtered.txt',
        expand('results/04_braken/{sample}.GTDB.centrifuge.k2report.T1.bracken.genus.report', sample=FID),
        expand('results/04_braken/{sample}.GTDB.centrifuge.k2report.T1.bracken.species.report', sample=FID),
        #expand('results/03_humann3Uniref50EC/{sample}_pathcoverage.tsv', sample=FID),
        'results/centrifuge.counts.all.txt',
        'results/centrifuge.counts.bracken.T1.genus.txt',
        'results/centrifuge.counts.bracken.T1.species.txt',



rule seqkitQCRaw:
    output:
        'results/00_qc/seqkit.report.raw.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a results/01_cutadapt/*.fastq.gz > {output} '



rule seqkitQCFiltered:
    output:
        'results/00_qc/seqkit.report.filtered.txt'
    conda:
        'seqkit'
    threads: 12
    shell:
        'seqkit stats -j {threads} -a results/02_kneaddata/*kneaddata.trimmed.fastq > {output} '



#TODO: consider using an input function to filter the seqkit summary tables for input to the profile.smk pipeline



localrules: generateCentrifugeSampleSheet



rule generateCentrifugeSampleSheet:
    output:
        sampleSheet='resources/centrifugeSampleSheet.tsv',
    threads:2
    shell: 
        './workflow/scripts/generate_centrifuge_sample_sheet.sh -d results/02_kneaddata -p kneaddata.trimmed.fastq -o {output.sampleSheet} '



rule centrifugeGTDB:
    input:
        sampleSheet='resources/centrifugeSampleSheet.tsv'
    output:
        out=expand('results/03_centrifuge/{sample}.GTDB.centrifuge', sample=FID),
        report=expand('results/03_centrifuge/{sample}.GTDB.centrifuge.report', sample=FID),
    log:
        'logs/centrifuge.GTDB.multi.log'
    conda:
        'centrifuge'
    threads: 32
    resources:
        mem_gb = lambda wildacards, attempt: 140 + ((attempt -1) + 20),
        time = "06:00:00"
    shell:
        'centrifuge ' 
        '-x /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/centrifuge/GTDB '
        '--sample-sheet {input.sampleSheet} '
        '-t '
        '--threads {threads} '
        '&> {log} '



rule centrifugeKrakenReport:
    input:
        centrifuge='results/03_centrifuge/{samples}.GTDB.centrifuge',
    output:
        centrifugeKraken2='results/03_centrifuge/{samples}.GTDB.centrifuge.k2report'
    log:
        'logs/{samples}.centrifuge.to.kraken2.log'
    conda:
        'centrifuge'
    threads: 2
    shell:
        'centrifuge-kreport '
        '-x /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/centrifuge/GTDB '
        '{input.centrifuge} > '
        '{output.centrifugeKraken2}'



rule brackenCentrifugeGenus:
    input:
        centrifugeKraken2='results/03_centrifuge/{samples}.GTDB.centrifuge.k2report',
    output:
        braken='results/04_braken/{samples}.GTDB.centrifuge.k2report.T1.bracken.genus',
        brakenReport='results/04_braken/{samples}.GTDB.centrifuge.k2report.T1.bracken.genus.report',
    log:
        'logs/{samples}.centrifuge.bracken.genus.GTDB.log'
    conda:
        'kraken2'
    threads: 2 
    shell:
        'bracken '
        '-d /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/kraken/GTDB '
        '-i {input.centrifugeKraken2} '
        '-o {output.braken} '
        '-w {output.brakenReport} '
        '-r 80 '
        '-l G '
        '-t 1 '
        '&> {log} '



rule brackenCentrifugeSpecies:
    input:
        centrifugeKraken2='results/03_centrifuge/{samples}.GTDB.centrifuge.k2report',
    output:
        braken='results/04_braken/{samples}.GTDB.centrifuge.k2report.T1.bracken.species',
        brakenReport='results/04_braken/{samples}.GTDB.centrifuge.k2report.T1.bracken.species.report',
    log:
        'logs/{samples}.centrifuge.bracken.species.GTDB.log'
    conda:
        'kraken2'
    threads: 2 
    shell:
        'bracken '
        '-d /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/kraken/GTDB '
        '-i {input.centrifugeKraken2} '
        '-o {output.braken} '
        '-w {output.brakenReport} '
        '-r 80 '
        '-l S '
        '-t 1 '
        '&> {log} '



rule humann3Uniref50EC:
    input:
        kneaddataReads='results/02_kneaddata/{samples}_kneaddata.trimmed.fastq'
    output:
        genes = 'results/03_humann3Uniref50EC/{samples}_genefamilies.tsv',
        pathways = 'results/03_humann3Uniref50EC/{samples}_pathabundance.tsv',
        pathwaysCoverage = 'results/03_humann3Uniref50EC/{samples}_pathcoverage.tsv'
    log:
        'logs/{samples}.human3.uniref50EC.log'
    conda:
        'biobakery'
    threads: 16
    resources:
        mem_gb= lambda wildcards, attempt: 24 + ((attempt - 1) + 12) 
    message:
        'humann3 profiling with uniref50EC: {wildcards.samples}\n'
    shell:
        'humann3 '
        '--memory-use minimum '
        '--threads {threads} '
        '--bypass-nucleotide-search '
        '--search-mode uniref50 '
        '--protein-database /bifo/scratch/2022-BJP-GTDB/biobakery/humann3/unirefECFilt '
        '--input-format fastq '
        '--output results/03_humann3Uniref50EC '
        '--input {input.kneaddataReads} '
        '--output-basename {wildcards.samples} '
        '--o-log {log} '
        '--remove-temp-output '



rule combineCentrifugeReports:
    input:
        expand('results/03_centrifuge/{sample}.GTDB.centrifuge.k2report', sample=FID),
    output:
        'results/03_centrifuge/rough.counts.all.txt'
    conda:
        'kraken2'
    threads: 2
    shell:
        'combine_kreports.py -o {output} -r {input} '



rule combineBrackenGenusReports:
    input:
        expand('results/04_braken/{sample}.GTDB.centrifuge.k2report.T1.bracken.genus', sample=FID),
    output:
        'results/centrifuge.counts.bracken.T1.genus.txt'
    conda:
        'kraken2'
    threads: 2
    shell:
        'combine_bracken_outputs.py -o {output} --files {input} '



rule combineBrackenSpeciesReports:
    input:
        expand('results/04_braken/{sample}.GTDB.centrifuge.k2report.T1.bracken.species', sample=FID),
    output:
        'results/centrifuge.counts.bracken.T1.species.txt'
    conda:
        'kraken2'
    threads: 2
    shell:
        'combine_bracken_outputs.py -o {output} --files {input} '



rule formatCombinedCentrifugeReport:
    input:
        'results/03_centrifuge/rough.counts.all.txt'
    output:
        'results/centrifuge.counts.all.txt'
    threads: 2
    shell:
        'workflow/scripts/reformat_centrifuge_count_matrix.sh -i {input} -p results/03_centrifuge && '
        'mv results/03_centrifuge/clean.count.matrix.txt {output} '
