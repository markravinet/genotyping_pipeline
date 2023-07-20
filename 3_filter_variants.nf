#!/usr/bin/env nextflow

// nextflow pipeline for variant filtering
// developed by Mark Ravinet - 07/06/2023
// v 1.2 - 20/06/2023

// script paramaters
anno_vcfs = Channel.fromPath( './vcf/*.vcf.gz' )

// Set filtering params
params.miss=0.8
params.q_site1=30
params.q_site2=30
params.min_depth=5
params.max_depth=30
params.min_geno_depth=5
params.max_geno_depth=30
params.keep="./path/to/file"

// Step 1 - normalise -  remove spanning deletions,indels and normalise
process normalise {

  input:
  path (anno_vcf)

  output:
  tuple \
    file ("${anno_vcf.simpleName}_norm.vcf.gz"), \
    file ("${anno_vcf.simpleName}_norm.vcf.gz.csi") 

  //path 'norm.vcf.gz*' into norm_vcfs

  """
  bcftools view -V indels -e 'ALT="*" | N_ALT>1' $anno_vcf | bcftools norm -D -O z -o ${anno_vcf.simpleName}_norm.vcf.gz

  bcftools index ${anno_vcf.simpleName}_norm.vcf.gz
  """

}

// Step 2 - filter for pop structure and genome scan
process filter_vcf {

  // publish simlinks into a final vcf directory
  publishDir 'vcf_filtered', saveAs: { filename -> "$filename" }, mode: 'copy'

  input:
  tuple file(norm_vcf), file(norm_vcf_index) 

  output:
  tuple \
    file ("${norm_vcf.simpleName}_filtered_ps.vcf.gz"), \
    file ("${norm_vcf.simpleName}_filtered_ps.vcf.gz.csi"), \
    file ("${norm_vcf.simpleName}_filtered_gs.vcf.gz"), \
    file ("${norm_vcf.simpleName}_filtered_gs.vcf.gz.csi") 

  """
  if [[ -f ${params.keep} ]]; then

    echo "File of individuals to filter provided - adding --keep option."

    # for pop structure
    vcftools --gzvcf $norm_vcf  --remove-indels --remove-filtered-all \
    --keep ${params.keep} \
    --min-alleles 2 --max-alleles 2 \
    --max-missing ${params.miss} --minQ ${params.q_site1} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout | \
    bcftools view -e 'N_ALT>1' -O z -o ${norm_vcf.simpleName}_filtered_ps.vcf.gz

    bcftools index ${norm_vcf.simpleName}_filtered_ps.vcf.gz

    # for genome scans
    vcftools --gzvcf $norm_vcf --remove-indels --remove-filtered-all \
    --keep ${params.keep} \
    --max-alleles 2 \
    --minQ ${params.q_site2} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout | \
    bcftools view -e 'N_ALT>1' -O z -o ${norm_vcf.simpleName}_filtered_gs.vcf.gz
    
    bcftools index ${norm_vcf.simpleName}_filtered_gs.vcf.gz

  else

    echo "Not filtering for specific individuals"
    # for pop structure
    vcftools --gzvcf $norm_vcf  --remove-indels --remove-filtered-all \
    --min-alleles 2 --max-alleles 2 \
    --max-missing ${params.miss} --minQ ${params.q_site1} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout | \
    bcftools view -e 'N_ALT>1' -O z -o ${norm_vcf.simpleName}_filtered_ps.vcf.gz

    bcftools index ${norm_vcf.simpleName}_filtered_ps.vcf.gz

    # for genome scans
      vcftools --gzvcf $norm_vcf --remove-indels --remove-filtered-all \
    --max-alleles 2 \
    --minQ ${params.q_site2} \
    --min-meanDP ${params.min_depth} --max-meanDP ${params.max_depth} \
    --minDP ${params.min_geno_depth} --maxDP ${params.max_geno_depth} \
    --recode --recode-INFO-all --stdout | \
    bcftools view -e 'N_ALT>1' -O z -o ${norm_vcf.simpleName}_filtered_gs.vcf.gz
    
    bcftools index ${norm_vcf.simpleName}_filtered_gs.vcf.gz

  fi
  """

}

// workflow starts here!

workflow{
  normalise(anno_vcfs) | filter_vcf
}