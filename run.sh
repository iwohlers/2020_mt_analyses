#!/bin/bash

echo "GENERATING AND ACTIVATING BIOCONDA WORKFLOW ENVIRONMENT..."
export PATH="/scratch/programs/miniconda3/bin:$PATH"
if [ ! -d /scratch/programs/miniconda3/envs/workflow_2020_population_genetics ]; then
    conda env create -n workflow_2020_population_genetics --file environment.yaml
fi
source activate workflow_2020_population_genetics

echo "RUNNING SNAKEMAKE WORKFLOW..."

snakemake -j 2 --rerun-incomplete -k -p --use-conda combined_mt_fasta combined_mt_vcf #combined_egypt_haplogroup_file combined_egypt_contamination_file fastqc_egyptian_summary get_coverage_egypt_all summarize_mapping_stats_egypt
