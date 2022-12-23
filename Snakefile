# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"

import os

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

print("Found: ")
for entry in FIDs:
    print(entry)



rule all:
    input:
        expand('04_kraken2GTDB/{samples}.GTDB.report.k2', samples = FIDs),
        # expand('03_kmcpGTDB/{samples}.search.tsv.gz', samples = FIDs),
        # expand('04_kraken2/{samples}.report.k2', samples = FIDs),

        # expand('03_humann/{samples}_kneaddata_pathabundance.tsv', samples = FIDs),
        # expand('03_kmcpGTDB/{samples}.profile.tsv', samples = FIDs),
        
        '00_qc/ReadsMultiQCReport.html',
        '00_qc/KDRReadsMultiQCReport.html'



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
        'Spewing barcodes for: {params.libraries}...\n'
    shell:
        'gquery '
        '-t gbs_keyfile '
        '-b library -p "columns=factid,barcode;fasta;noheading;no_unpivot" '
        '{params.libraries} > '
        '{output.barcodes} 2> '
        '{log}'



rule cutadapt: #demultiplexing GBS reads
    output:
        expand('01_cutadapt/{samples}.fastq.gz', samples = FIDs),
    input:
        barcodes = rules.generateBarcodes.output.barcodes,
        lane01 = config['novaseq']['lane01'],
        lane02 = config['novaseq']['lane02'],
    conda:
        'cutadapt'
        # 'docker://quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'
    threads: 16
    resources:
        mem_gb=32,
    message:
        'Demultiplexing lanes...'
    shell:
        'zcat {input.lane01} {input.lane02} | '
        'cutadapt '
        '-j 16 '
        '--discard-untrimmed '
        '--no-indels '
        '-g ^file:{input.barcodes} '
        r'-o "01_cutadapt/{{name}}.fastq.gz" '
        '-' # indicates stdin



rule fastqc:
    output:
        html = '00_qc/fastqc/{samples}_fastqc.html',
        zip = "00_qc/fastqc/{samples}_fastqc.zip"
    input:
        fastq = '01_cutadapt/{samples}.fastq.gz'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
    message:
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o 00_qc/fastqc/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'



rule multiQC:
    output:
        multiQC='00_qc/ReadsMultiQCReport.html'
    input:
        fastqc= expand('00_qc/fastqc/{samples}_fastqc.zip', samples = FIDs)
    conda:
        'multiqc'
        # 'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '-n 00_qc/ReadsMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'



#TODO Rule to Build Rambv2 index



rule kneaddata:
    output:
        trimReads = temp('02_kneaddata/{samples}_kneaddata.trimmed.fastq'),
        trfReads = temp('02_kneaddata/{samples}_kneaddata.repeats.removed.fastq'),
        ovineReads = temp('02_kneaddata/{samples}_kneaddata_GCF_016772045.1-ARS-UI-Ramb-v2.0_bowtie2_contam.fastq'),
        KDRs = '02_kneaddata/{samples}_kneaddata.fastq',
        readStats = '02_kneaddata/{samples}.read.stats.txt'
    input:
        reads = '01_cutadapt/{samples}.fastq.gz',
    conda:
        'biobakery'
    log:
        'logs/kneaddata/{samples}.kneaddata.log'
    threads: 4
    resources:
        mem_gb=8,
        time='02:00:00'
    message:
        'kneaddata: {wildcards.samples}\n'
    shell:
        'kneaddata '
        '--trimmomatic-options "MINLEN:60 ILLUMINACLIP:/home/perrybe/conda-envs/biobakery/share/trimmomatic-0.39-2/adapters/illuminaAdapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 CROP:80" '
        '--input {input.reads} '
        '-t {threads} '
        '--log-level INFO '
        '--log {log} '
        '--trimmomatic /home/perrybe/conda-envs/biobakery/share/trimmomatic '
        '--sequencer-source TruSeq3 '
        '-db ref/Rambv2/GCF_016772045.1-ARS-UI-Ramb-v2.0 '
        '-o 02_kneaddata && '
        'seqkit stats -j 12 -a 02_kneaddata/{wildcards.samples}*.fastq > {output.readStats}'



#TODO Compress output reads



rule fastqcKDRs:
    output:
        "00_qc/fastqcKDR/{samples}_kneaddata_fastqc.zip"
    input:
        fastq = '02_kneaddata/{samples}_kneaddata.fastq'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
    message:
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o 00_qc/fastqcKDR/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'



rule multiQCKDRs:
    output:
        '00_qc/KDRReadsMultiQCReport.html'
    input:
        fastqc= expand('00_qc/fastqcKDR/{samples}_kneaddata_fastqc.zip', samples = FIDs)
    conda:
        'multiqc'
        # 'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '-n 00_qc/KDRReadsMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'



rule vsearchUniques:
    output:
        uniqueReads="03_uniques/{samples}.uniques.merged.fastq.gz",
    input:
        KDRs=rules.kneaddata.output.KDRs,
    log:
        "logs/vsearchUniques/{samples}.vsearch.unqiues.log",
    conda:
        "vsearch"
    threads: 1
    message:
        "dereplicating: {wildcards.samples}\n"
    shell:
        "vsearch "
        "--gzip "
        "--threads {threads} "
        "--log {log} "
        "--fastx_uniques {input.KDRs} "
        "--sizein "
        "--minuniquesize 1 "
        "--relabel_self "
        "--sizeout "
        "--fastqout - | gzip > {output.uniqueReads}"



rule kraken2:
    output:
        k2Out='04_kraken2/{samples}.out.k2',
        k2Report='04_kraken2/{samples}.report.k2'
    input:
        KDRs=rules.vsearchUniques.output.uniqueReads
    log:
        'logs/{samples}.kraken2.log'
    conda:
        'kraken2'
    threads: 10 
    resources: 
        mem_gb=146,
        partition="inv-bigmem,inv-bigmem-fast"
    shell:
        'kraken2 '
        '--db ref/kraken2 '
        '--report {output.k2Report} '
        '--report-minimizer-data '
        '{input.KDRs} > {output.k2Out}'


rule GTDBtoRam:
    input:
        GTDBdir='/bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/kraken/GTDB',        
        # k2hash='/bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/kraken/GTDB/hash.k2d',
        # k2opts='/bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/kraken/GTDB/opts.k2d',
        # k2taxo='/bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/kraken/GTDB/taxo.k2d',
    output:
        kraken2GTDB=temp(directory('/dev/shm/GTDB')),
        # ramHash=temp('/dev/shm/GTDB/hash.k2d'),
        # ramOpts=temp('/dev/shm/GTDB/opts.k2d'),
        # ramTaxo=temp('/dev/shm/GTDB/taxo.k2d'),
    conda:
        'kraken2'
    threads: 4
    resources:
        partition="inv-bigmem"
    shell:
        # 'mkdir -p {output.kraken2GTDB}; '
        # 'cp {input.k2hash} {output.ramHash}; '
        # 'cp {input.k2opts} {output.ramOpts}; '
        # 'cp {input.k2taxo} {output.ramTaxo}; '
        'echo "starting copy: $(date)... '
        'cp -p {input.GTDB}/*.k2d {output.kraken2GTDB} ; '
        'echo "copy completed: $(date)... '
        'ls -lh /dev/shm/GTDB; '
        'du -sh /dev/shm/GTDB; '



rule kraken2GTDB:
    output:
        k2OutGTDB='04_kraken2GTDB/{samples}.GTDB.out.k2',
        k2ReportGTDB='04_kraken2GTDB/{samples}.GTDB.report.k2'
    input:
        KDRs=rules.vsearchUniques.output.uniqueReads,
        kraken2GTDB=rules.GTDBtoRam.output.kraken2GTDB
    log:
        'logs/{samples}.kraken2.GTDB.log'
    conda:
        'kraken2'
    threads: 4 
    resources: 
        partition="inv-bigmem"
    shell:
        'kraken2 '
        '--db {input.kraken2GTDB} '
        '--memory-mapping '
        '--report {output.k2ReportGTDB} '
        '--report-minimizer-data '
        '{input.KDRs} > {output.k2OutGTDB}'



# rule braken:
#     output:
#         '04_kraken2/{samples}.braken.out'
#     input:
#         k2out=rules.kraken2.output
#     log:
#         'logs/{samples}.braken.log'
#     conda:
#         ''# TBD
#     threads: 8
#     message:
#         ''
#     shell:
#         ''



rule humann3:
    output:
        genes = '03_humann/{samples}_kneaddata_genefamilies.tsv',
        pathways = '03_humann/{samples}_kneaddata_pathabundance.tsv',
        pathwaysCoverage = '03_humann/{samples}_kneaddata_pathcoverage.tsv'
    input:
        KDRs=rules.kneaddata.output.KDRs
    log:
        'logs/{samples}.human3.log'
    conda:
        'biobakery'
    threads: 8
    resources:
        mem_gb=24,
    message:
        'humann3 profiling: {wildcards.samples}\n'
    shell:
        'humann '
        '--threads {threads} '
        '--input {input.KDRs} '
        '--output 03_humann '
        '--nucleotide-database ref/biobakery/humann_dbs/chocophlan '
        '--bypass-prescreen '
        #'--bypass-nucleotide-search '
        '--memory-use maximum '
        '--input-format fastq '
        '--search-mode uniref90 '
        '--protein-database ref/biobakery/humann_dbs/uniref ' #TODO Update
        '--verbose '
        '--log-level INFO '
        '--o-log {log} '
        '--remove-temp-output '


rule kmcpSearch:
    output:
        search='03_kmcpGTDB/{samples}.search.tsv.gz'
    input:
        kmcpGTDB='/dataset/2022-BJP-GTDB/active/kmcp/gtdb.kmcp',
        KDRs=rules.vsearchUniques.output.uniqueReads,
    log:
        'logs/{samples}.kmcp.search.GTDB.log'
    conda:
        'kmcp'
    threads:12
    resources: 
        mem_gb=104,
        partition="inv-bigmem,inv-bigmem-fast,inv-iranui-fast,inv-iranui"
    shell:
        'kmcp search '
        '-w '
        '--threads {threads} '
        '--db-dir {input.kmcpGTDB} '
        '{input.KDRs} '
        '-o {output.search} '
        '--log {log}'     



rule kmcpProfile:
    output:
        profile='03_kmcpGTDB/{samples}.profile.tsv'
    input:
        kmcpSearch=rules.kmcpSearch.output.search,
        kmcpTaxdump='/dataset/2022-BJP-GTDB/active/kmcp/gtdb-taxdump/R207',
        taxid='/dataset/2022-BJP-GTDB/active/kmcp/taxid.map',
    log:
        'logs/{samples}.kmcp.profile.log'
    threads: 8 
    resources: 
        mem_gb=6,
    shell:
        'kmcp profile '
        '--mode 1 '
        '--threads {threads} '
        '-X {input.kmcpTaxdump} '
        '-T {input.taxid} '
        '-o {output.profile} '
        '--log {log} '
        '{input.kmcpSearch} '



# # rule human3GTDB:
#     output:

#     input:

#     log:

#     threads:

#     message:

#     shell:

