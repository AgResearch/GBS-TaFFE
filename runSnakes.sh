#!/bin/bash
### Wrapper script for launching snakemake workflows ###
echo "Launching snakemake workflow..."
sleep 3
echo "Beginning executing on: $(date)"

source activate snakemake

echo "Demultiplexing SQ libraries"
#snakemake --profile config/slurm --snakefile workflow/demux.smk
echo "Demultiplexing completed"

echo "Preparing fastq for profiling."

echo "Collecting SQ runs: $(ls results/SQ*).."


snakemake --profile config/slurm --snakefile workflow/prepare.smk
echo "Preparations completed."

echo "Starting profiling..."
snakemake --profile config/slurm --snakefile workflow/profile.smk
echo "profiling completed."

##snakemake --rulegraph | dot -T svg > rulegraph.svg

conda deactivate

echo "snakemake run completed: $(date)"
