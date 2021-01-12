#!/bin/bash

echo "GENERATING AND ACTIVATING BIOCONDA WORKFLOW ENVIRONMENT..."
export PATH="/scratch/programs/miniconda3/bin:$PATH"
if [ ! -d /scratch/programs/miniconda3/envs/workflow_2020_population_genetics ]; then
    conda env create -n workflow_2020_population_genetics --file environment.yaml
fi
source activate workflow_2020_population_genetics

echo "RUNNING SNAKEMAKE WORKFLOW..."

snakemake --rerun-incomplete -p --use-conda combined_egypt_haplogroup_file combined_egypt_contamination_file
