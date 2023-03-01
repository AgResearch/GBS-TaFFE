#!/bin/bash

### Wrapper script for launching snakemake workflows ###
clear
echo "Launching snakemake workflow..."
sleep 3
echo "Beginning executing on: $(date)"

source activate snakemake

snakemake --profile config/slurm

snakemake --report snakemakeReport.html

snakemake --rulegraph | dot -T svg > rulegraph.svg

source deactivate

echo "snakemake run completed: $(date)
