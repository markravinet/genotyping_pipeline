#!/bin/bash

# Job name:
#SBATCH --job-name=concat_vcf

# partition/queue job being run on
#SBATCH --partition=defq

# number of nodes
#SBATCH --nodes=1

# tasks per node
#SBATCH --ntasks-per-node=1

# cpus-per-task
#SBATCH --cpus-per-task=8

# memory - max memory job can use
#SBATCH --mem=12g

# total time - i.e. wallclock limit - how long job has to run
#SBATCH --time=24:00:00

# notify job failure
#SBATCH --mail-user=mark.ravinet@nottingham.ac.uk
#SBATCH --mail-type=FAIL

## set up job environment

# define variables - nb TARGET is defined from command line
VCF_FILE=/path/to/list/of/vcfs
OUTPUT_VCF1=/nextflow_workdir/vcf/sparrows_variants.vcf.gz
OUTPUT_VCF2=/nextflow_workdir/vcf/sparrows_variants_norm.vcf.gz

# concat vcfs
bcftools concat -f $VCF_FILE --threads 8 -n -O z -o $OUTPUT_VCF1

# index the sorted vcf1
bcftools index $OUTPUT_VCF1

# next normalise
bcftools norm -d none -O z -o $OUTPUT_VCF2 $OUTPUT_VCF1

bcftools index $OUTPUT_VCF2
