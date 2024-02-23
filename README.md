# Analysing RE-RRSeq with GTDB taxonomic framework


# Pipeline Overview

## Demultiplexing GBS libraries for analysis with GBS-TaFFE
Prior to running preprocessing and classification using the maink GBS-TaFFE.smk workflow, the library level fastq file(s) must first be demultiplexed. Library fastq files have the convention SQ#### as an identifier. Currently the demux.smk workflow is used to run cutadapt to demultiplex a library. This depends on two inputs: 1) a gbs cutadapt formatted barcodes file generated using gquery, 2) a symbolic link of the library in the /fastq directory (often this is a symlink to the link farm; the link farm is populated via the gbs_prism pipeline).

```
GBS-TaFFE ~ $ tree fastq/ | head
fastq/
├── README.md
├── SQ0738_CCHC7ANXX_s_1_fastq.txt.gz -> /dataset/hiseq/active/fastq-link-farm/SQ0738_CCHC7ANXX_s_1_fastq.txt.gz
├── SQ0738_CCHC7ANXX_s_2_fastq.txt.gz -> /dataset/hiseq/active/fastq-link-farm/SQ0738_CCHC7ANXX_s_2_fastq.txt.gz
├── SQ0739_CCHC7ANXX_s_3_fastq.txt.gz -> /dataset/hiseq/active/fastq-link-farm/SQ0739_CCHC7ANXX_s_3_fastq.txt.gz
├── SQ0739_CCHC7ANXX_s_4_fastq.txt.gz -> /dataset/hiseq/active/fastq-link-farm/SQ0739_CCHC7ANXX_s_4_fastq.txt.gz
├── SQ0740_CCHC7ANXX_s_5_fastq.txt.gz -> /dataset/hiseq/active/fastq-link-farm/SQ0740_CCHC7ANXX_s_5_fastq.txt.gz
├── SQ0740_CCHC7ANXX_s_6_fastq.txt.gz -> /dataset/hiseq/active/fastq-link-farm/SQ0740_CCHC7ANXX_s_6_fastq.txt.gz
├── SQ0741_CCHC7ANXX_s_7_fastq.txt.gz -> /dataset/hiseq/active/fastq-link-farm/SQ0741_CCHC7ANXX_s_7_fastq.txt.gz
├── SQ0741_CCHC7ANXX_s_8_fastq.txt.gz -> /dataset/hiseq/active/fastq-link-farm/SQ0741_CCHC7ANXX_s_8_fastq.txt.gz
...

GBS-TaFFE ~ $ tree resources/ | head
resources/
├── SQ0738.cutadapt.barcodes.fasta
├── SQ0739.cutadapt.barcodes.fasta
├── SQ0740.cutadapt.barcodes.fasta
├── SQ0741.cutadapt.barcodes.fasta
...
```

