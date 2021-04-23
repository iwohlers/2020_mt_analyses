#!/bin/bash

echo "GENERATING AND ACTIVATING BIOCONDA WORKFLOW ENVIRONMENT..."
export PATH="/scratch/programs/miniconda3/bin:$PATH"
if [ ! -d /scratch/programs/miniconda3/envs/workflow_2020_population_genetics ]; then
    conda env create -n workflow_2020_population_genetics --file environment.yaml
fi
source activate workflow_2020_population_genetics

echo "RUNNING SNAKEMAKE WORKFLOW..."

snakemake -j 8 --rerun-incomplete -k -p --use-conda  combined_mt_fasta_1000g combined_mt_fasta_hgdp combined_mt_vcf_1000g combined_mt_vcf_hgdp #extract_extra_high_quality_incl_minor  #north_african_mt/mitoimpute_plus_north_african_mt.fa #combined_sudan_haplogroup_file combined_sudan_contamination_file fastqc_summary get_coverage_all summarize_mapping_stats
