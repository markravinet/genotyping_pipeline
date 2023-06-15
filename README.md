# Variant calling pipeline

### Mark Ravinet 
### 13/06/2023

## Introduction

Welcome to the readme and guide for the Sparrow Ecological & Evolutionary Genommics group variant calling pipeline. After many years of calling variants using a set of multiple slurm scripts, I decided to write a nextflow pipeline that vastly simplifies the process, meaning that you will only have to submit three scripts to go from raw reads to a final, filtered `vcf` ready for analysis.

The scripts are designed to require minimal input from a user but they will produce finished output, statistics and reports for what programs where used in the process. As well as simplifying the process, these scripts are designed to **standardise** the variant calling used within the group, i.e. to prevent different projects using different pipelines and leading to incompatible datasets. However, more importantly, the emphasis here is to improve **reproducibility** - i.e. promoting open science and making it possible for anyone interested to reproduce results in a straightforward way.

There are three steps in the pipeline, each run using a specific nextflow script:

1. `trim_map_realign` - trim reads for quality and adapter sequences, map to a reference genome and perform indel realignment
2. `call_variants` - call variants across genome windows, concatenate them together to produce per chromosome vcf files
3. `filter_variants` - filter vcfs for population structure and genome scan analysis

## Installation

The simplest way to install everything you need for these scripts to work is to use template `conda` environment. This means you first need to install and setup `conda`. 

### Installing conda

The full `conda` installation is not necessary. Instead you can use miniconda. Go here and copy the link address for the latest release for Linux 64-bit. Then use `wget` to download it like so:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

From your home directory, simply run this script:

```
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the prompts and this will install `conda`. Once you have done that, you will need to update it to ensure there are no issues.

```
conda update conda
```

We are then read to set up the `conda` environment you need to run the scripts.

### Configuring the conda environment



Finally in order to ensure `abra2` can access the libraries it needs to run, you will need to add the following line to your `bash_profile`:

```
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/lib
```

## Step 1  - trimming, mapping and aligning

This first script will take your raw reads, trim them for low-quality bases and remove any adapter sequences. It will then map them to a reference genome of your choice (default is the 2014 House sparrow reference). Following this, it will perform realignment around indels using `abra2` to ensure that short insertions or deletions do not cause false positive variant calls. It also ensures that the samples are renamed to the correct name for all downstream analysis. Finally, the script will produce statistics on the mapping efficiency and depth of coverage of each mapped individual. 

### The input csv format

In order to run this script, you need to provide it with a csv file with the following columns:

1. Original sample name (only necessary for older samples, otherwise this can just be the sample name)
2. New sample name (i.e. UoN code) - if this is already the sample name in 1, then just repeat it here
3. Forward read location - this should be the **full path** to the forward read
4. Reverse read location - this should be the **full path** to the reverse read
5. Adapter name - the name for the adapter sequences - this should be one of the set of options outlined below.

Adapters should be written exactly as one of the following options:
- `Illumina_UD-PE`
- `NexteraPE-PE`
- `TruSeq2-PE`
- `TruSeq3-PE-2`
- `TruSeq3-PE`
- `TruSeqUD-PE`

**Note** If your individual has multiple forward and reverse reads, these should all be entered into a separate row for the csv (i.e. for each forward and reverse pair). This will happen when an individual is sequenced across multiple lanes for example. The script is built to account for this and you can give each line the same sample name - **however** the read locations must be different for each. 

An example of the file format is shown below:

| old_sample1_name | new_sample1_name | /path/to/forward_read1/ | /path/to/reverse_read1/ | TruSeq3-PE |
| old_sample2_name | new_sample2_name | /path/to/forward_read2/ | /path/to/reverse_read2/ | TruSeqUD-PE |

### Running the script

Once you have your input csv file ready, you can run the script. To do this, you just need to do the following:

```
nextflow run 1_trim_map_realign_v0.2.nf --samples samples_test.csv
```

Note that if you run the script in this way, it will run in default mode and use the house sparrow reference genome. In order to alter that, you can add an additional option, `--ref` and specify the location of this alternative reference. Note that in order to do that, you **must** ensure that reference genome has been indexed by bwa prior to running the analysis. There is info on [how to do that here](https://speciationgenomics.github.io/mapping_reference/):

```
nextflow run 1_trim_map_realign_v0.2.nf --samples samples_test.csv --ref /path/to/alt/ref_genome.fa
```
### Resuming scripts if there are issues

On occassion your script might run into issues. For tips on how to resolve these, see the **Troubleshooting** section below. Once you have fixed this, you can actually restart the nextflow script from a specific point using the `-resume` command. For example:

```
nextflow run 1_trim_map_realign_v0.2.nf --samples samples_test.csv -resume
```

Remember the script might not pick up exactly where it left off as it runs from preset checkpoints. However, it will generally be much faster doing this than starting from scratch.

### Script outputs

Once complete, the script will create two directories with outputs in:

1. `align` - this contains the mapped, realigned and sorted bamfiles with their indexes for each individual you ran the script on.
2. `stats` - this directory contains statistics from each of the bamfiles for their mapping success and depth of coverage. More information below.



## Step 2  - trimming, mapping and aligning

## Step 3