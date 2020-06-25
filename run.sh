#!/bin/bash

echo "GENERATING AND ACTIVATING BIOCONDA WORKFLOW ENVIRONMENT..."
export PATH="/scratch/programs/miniconda3/bin:$PATH"
if [ ! -d /scratch/programs/miniconda3/envs/workflow_2020_mt_analyses ]; then
    conda env create -n workflow_2020_mt_analyses --file environment.yaml
fi
source activate workflow_2020_mt_analyses

echo "RUNNING SNAKEMAKE WORKFLOW..."

snakemake --rerun-incomplete -p -j 12 --use-conda combined_haplogroup_file
