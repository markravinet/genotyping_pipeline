# Variant calling pipeline

### Mark Ravinet 
### 20/06/2023

## Introduction

Welcome to the readme and guide for the Sparrow Ecological & Evolutionary Genommics group variant calling pipeline. After many years of calling variants using a set of multiple slurm scripts, I decided to write a nextflow pipeline that vastly simplifies the process, meaning that you will only have to submit three scripts to go from raw reads to a final, filtered `vcf` ready for analysis.

The scripts are designed to require minimal input from a user but they will produce finished output, statistics and reports for what programs where used in the process. As well as simplifying the process, these scripts are designed to **standardise** the variant calling used within the group, i.e. to prevent different projects using different pipelines and leading to incompatible datasets. However, more importantly, the emphasis here is to improve **reproducibility** - i.e. promoting open science and making it possible for anyone interested to reproduce results in a straightforward way.

There are three steps in the pipeline, each run using a specific nextflow script:

1. `trim_map_realign` - trim reads for quality and adapter sequences, map to a reference genome and perform indel realignment.
2. `call_variants` - call variants across genome windows, concatenate them together to produce per chromosome vcf files.
3. `filter_variants` - filter vcfs for population structure and genome scan analysis.

In addition to these scripts, there are three helper scripts which are optional. These are:

1. `create_genome_windows` - a shell script to split your genome of choice into windows for `call_variants` to run on. So you can split the genome into 10 Mb windows and call variants much more efficiently.
2. `concat_vcfs` - concatenate the variant vcfs created by `filter_variants` into a single vcf for the entire genome

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

An example of the file format is shown below - note that the true file **should not have headers**:
| old name | new name | read1 path | read2 path | adapter |
|------------------|------------------|-------------------------|-------------------------|------------|
| old_sample1_name | new_sample1_name | /path/to/forward_read1/ | /path/to/reverse_read1/ | TruSeq3-PE |
| old_sample2_name | new_sample2_name | /path/to/forward_read2/ | /path/to/reverse_read2/ | TruSeqUD-PE |

### Running the script

Once you have your input csv file ready, you can run the script. To do this, you just need to do the following:

```
nextflow run 1_trim_map_realign.nf --samples samples_test.csv
```

Note that if you run the script in this way, it will run in default mode and use the house sparrow reference genome. In order to alter that, you can add an additional option, `--ref` and specify the location of this alternative reference. Note that in order to do that, you **must** ensure that reference genome has been indexed by bwa prior to running the analysis. There is info on [how to do that here](https://speciationgenomics.github.io/mapping_reference/):

```
nextflow run 1_trim_map_realign.nf --samples samples_test.csv --ref /path/to/alt/ref_genome.fa
```
### Resuming scripts if there are issues

On occassion your script might run into issues. For tips on how to resolve these, see the **Troubleshooting** section below. Once you have fixed this, you can actually restart the nextflow script from a specific point using the `-resume` command. For example:

```
nextflow run 1_trim_map_realign.nf --samples samples_test.csv -resume
```

Remember the script might not pick up exactly where it left off as it runs from preset checkpoints. However, it will generally be much faster doing this than starting from scratch.

### Script outputs

Once complete, the script will create two directories with outputs in:

1. `align` - this contains the mapped, realigned and sorted bamfiles with their indexes for each individual you ran the script on.
2. `stats` - this directory contains statistics from each of the bamfiles for their mapping success and depth of coverage. More information below.

## Step 2 - Variant calling

The second script in the pipeline will take a list of bamfiles and performs genotyping on all individuals. To do this, it uses `bcftools` and will call sites at every position in the genome (i.e. it calls invariant sites as well as variants). This is obviously a large job, especially on larger genomes. So to increase efficiency, the script parallelises across genome windows. The default is 10 Mb but you can set these to whatever size you wish. However, tweaking windows has to be done with a separate bash script (see below), not within the nextflow pipeline. After calling genotypes in windows, the script will take care of sorting and concatenating the windows together so that you are left with a vcf file for each chromosome, the mtDNA and also the unanchored scaffolds in your genome.

### Setting up genome windows

In order to parallelise across genome windows, the script requires a list of said windows in a text file. These are easy to generate using the helper script `0_create_genome_windows.sh`. This simple bash script uses bedtools to split the genome into windows of whatever size you wish (I recommend 10,000,000 bases for most bird genomes) and then it will also create separate files to account for all the scaffolds.

You need to run this script in the directory you are running the nextflow analysis in. You run it like so:

```
bash 0_create_genome_windows.sh ref_index window_size output_name
```

A more specific example, that will work for the house sparrow genome looks like this:

```
bash 0_create_genome_windows.sh /share/Passer/data/reference/house_sparrow_ref.fa.fai 10000000 sparrow_genome_windows.list
```

So option 1 is the path to the reference index, in `fai` format, option 2 is the window size (10 Mb here) and option 3 is the name of the output.

Running this script will produce a set of different files. The first will be text file with a list of all the windows, i.e. `sparrow_genome_windows.list` in the example above. The second will be 10 files called `scaffolds:00`, `scaffolds:01` and so on. The pipeline needs all these files to be present in the base directory you are running nextflow in. 

### Creating a list of bams

The other input this script needs is a list of bamfiles. This is very simple - it is just a list of paths of the files that you intend to analyse. If you generated these using `1_trim_map_realign.nf` then the names should already be standardised as `samplename_realigned.bam`. Provided your bams are named this way, then the calling script will also ensure that the sample names are written into the final vcf. Here is an example of what the file should look like:

```
/path/to/align/PDOMNOR8934547_realigned.bam
/path/to/align/PDOMNOR8L19766_realigned.bam
/path/to/align/PDOMNOR8L19786_realigned.bam
/path/to/align/PDOMNOR8L52141_realigned.bam
/path/to/align/PDOMNOR8L52830_realigned.bam
```

Note that these bams do not need top be in the same directory and they do not need to be in a directory called align. This means you can call a vcf from bams in multiple locations easily. 

### Running the script

Once you have your list of bams and genome windows file ready, you can run the script like this:

```
nextflow run 2_call_variants.nf --bams bams.list --windows sparrow_genome_windows.list
```

As with all the scripts, you can use the `-resume` option to rerun from a checkpoint if it fails for any reason.

### Script outputs

The outputs for this script are much simpler than the previous step - it will create a directory called `vcf` and inside will be the gz compressed vcf files for each chromosome, the mtDNA and the genome unanchored scaffolds. These will be raw (i.e. unfiltered) and will contain calls for all variant and non-variant sites in the genome. The directory will also include the indexes for these vcfs.

## Step 3 - filtering

The final script takes control of filtering your vcf files and prepares them for downstream analysis. First it normalises them to remove any issues from concatenating across windows in the calling step. Then it applies custom filters using `vcftools` to create two sets of vcfs; one for population structure analyses (i.e. variants only) and one for genome scan analyses (i.e. variant and invariant sites). There is more information on the outputs below.

### Running the script

The filtering script is the easiest of the three main pipeline scripts to run. It does not require any input as it will automatically look for any vcfs which are gzipped, indexed and which are stored in a directory called `./vcf` in the base directory it is run in. This means it can simply be run with the default filtering options like so:

```
nextflow run 3_filter_variants.nf 
```

However, [as shown here](https://speciationgenomics.github.io/filtering_vcfs/), it is not a good idea to just run filters without checking whether they apply to your dataset. Instead you are able to tweak the filters with a number of options. These are simply provided to the script using option flags and are modified versions of the options for [vcftools](https://vcftools.github.io/examples.html)

- `--miss` - set the missing data at a value between 0 and 1 (where 0 allows 100% missing data and 1 means no missing data); default is 0.8
- `--q_site1` - site quality threshold (as a phred score) for the population structure vcfs - default is 30
- `--q_site2` - site quality threshold (as a phred score) for the genome scan vcfs - default is 30
- `--min_depth` - minimum mean depth of coverage for a variant across all samples - default is 5
- `--max_depth` - maximum mean depth of coverage for a variant across all samples - default is 30
- `--min_geno_depth` - minimum genotype depth per sample. If lower than this value, the genotype will be converted to a missing site - default is 5
- `--max_geno_depth` - maximum genotype depth per sample. If lower than this value, the genotype will be converted to a missing site - default is 30

You can provide all or some of these options to the script using these options. A fully worked example is below:

```
nextflow run 3_filter_variants.nf --miss 0.5 --q_site1 30 --q_site2 40 --min_depth 5 --max_depth 15 --min_geno_depth 5 --max_geno_depth 15
```

You do not need to provide all the options - for example, if you want to just alter the missing data threshold, the following will work.

```
nextflow run 3_filter_variants.nf --miss 0.5 
```

### Script outputs

The script will create a directory called `vcf_filtered`. Inside this vcf will be per chromosome (and scaffold) vcfs with two different suffixes.

- `chrXX_norm_filtered_ps.vcf.gz` - this is the population structure analysis vcf - it contains variant biallelic SNPs only - ready for PCA, ADMIXTURE and so on.
- `chrXX_norm_filtered_gs.vcf.gz` - this is the genome scan vcf - it contains variant and invariant sites - it is ready for both selection and introgression scans

### Creating a variants whole genome vcf

The `3_filter_variants.nf` script will produce outputs per chromosome. This is done for downstream efficiency but it is sometimes useful to have a whole genome vcf for variants - i.e. the `chrXX_norm_filtered_ps.vcf.gz` vcfs. To produce this, you can use the script `4_concat_vcfs.slurm` which will concatenate all the chromosomes together and normalise the output. This is a slurm script that can be submitted but it needs editing first. You should edit the three variables declared at the start of the script:

- `VCF_FILE` - a text file with the paths of the vcfs to concatenate into a single vcf - this MUST be in order
- `OUTPUT_VCF1` - the name of the output vcf prior to normalisation
- `OUTPUT_VCF2` - the name of the output vcf following normalisation

The `VCF_FILE` should look like this:

```
/path/to/my/vcfs/chr1__norm_filtered_ps.vcf.gz
/path/to/my/vcfs/chr2__norm_filtered_ps.vcf.gz
/path/to/my/vcfs/chr3__norm_filtered_ps.vcf.gz
/path/to/my/vcfs/chr4__norm_filtered_ps.vcf.gz
/path/to/my/vcfs/chr5__norm_filtered_ps.vcf.gz
```

The order is important here and should match the order of the chromosomes in the reference genome. 

