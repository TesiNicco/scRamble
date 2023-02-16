# scRamble
scRamble is a pipeline written in R to scramble chromosomes across a set of individuals to avoid concerns about privacy of genetic data.
This should be done prior to submission to non-European imputation panels due to privacy laws in the European union.

# In a nutshell
Due to the recent GDPR regulations, it is not possible anymore to share genetic data outside European union. This may represent an issue for sharing data with collaborators outside Europe, or to run applications on servers relying outside Europe. One of the most common application when dealing with genetic data from SNP array is genotype imputation. Genotype imputation is hosted in a server in The United States (https://imputation.biodatacatalyst.nhlbi.nih.gov/). In order to use the server, users need to upload genetic data. This creates a conflict with the current GDPR regulation. To address this, we developed scRamble, a simple pipeline that shuffles chromosomes across individuals making it impossible to reconstruct an individual's genome without the right key. scRamble takes genotype data, creates random identifiers for each sample, and shuffles individual's data across chromosome, while keeping locally the key to reassemble the chromosomes for each individual. After scRamble has shuffled genotypes, the input can be used for genotype imputation. As genotypec imputation is done at chromosome level, the approach doesn't compromise quality of the imputed data. After data is imputed, scRamble reconstruct the original chromosomal order for each individual.

# Installation
Clone `scRamble` package into your system by typing `git clone https://github.com/TesiNicco/scRamble.git`.
In order to run, `scRamble` needs the following packages:
1. `R` with the following libraries installed: `data.table` and `stringr`. You can easily install these libraries by opening `R` and typing `install.packages('data.table')` and `install.packages('stringr')`
2. bcftools (https://samtools.github.io/bcftools/bcftools.html)
3. plink (version 1.9 and version 2, available at https://www.cog-genomics.org/plink/)

You can install all packages with `Conda`. Please use the `.yml` file provided to install all necessary packages. To do so, make sure you have conda available in your system (https://docs.conda.io/en/latest/), and type:
1. Create conda environment: `conda env create --name scRamble --file=scRamble/bin/scRamble_env.yml`
2. Activate conda environment: `conda activate scRamble`
3. You may still need to install `R` packages manually.

# Input data
Input data for scRamble is a single VCF file including all SNPs and all individuals. The output directory where to place scrambled genotypes should be also provided.

# How to run
To run scRamble, type:
`Rscript path/to/scRamble/bin/scRamble.R [path/to/input.vcf.gz] [path/to/output_directory_name] [consider sex chromosome (yes/no)] [reference genome (hg19/hg38)]`

For example, assuming that (i) the input file is called `chrAll.input.vcf.gz`, (ii) you want the output files to be placed in a folder `scrambled_genotypes`, you are interested in sex chromosomes and your genome version is hg38 then type:
`Rscript bin/scRamble_sexchr.R chrAll.input.vcf.gz scrambled_genotypes yes hg38`

It's important to note that `scRamble` assumes data to be aligned with the reference genome of interest, and will not perform liftover of the data.

Depending on the run parameters, it may take some minutes to perform.

# Output scrambled data
After the run, you will find in the desired folder three output files:
1. the scrambled file in VCF format, one per chromosome. Each chromosome-file is sorted and indexed.
2. mapping file to be used for unscrambling [`Mapping_file.txt`]

# Reconstruct original data after imputation
After genotype imputation is performed, you want to unscramble the outputs and reconstruct the original genotypes.
To do so, move into the scRamble directory and type:
`Rscript path/to/bin/unscRamble.R [path/to/Mapping_file.txt] [path/to/chr1_imputed.vcf.gz] [consider sex chromosome (yes/no)] [path/to/output_directory]`

For example, assuming that (i) mapping file is called Mapping_file.txt, (ii) the imputed file, usually already divided in chromosomes, 
is called chr1_imputed.vcf.gz (make sure to specificy chromosome 1 here, as the others will be infer from here), (iii) you want to consider sex chromosomes, and (iv) the output directory is output_directory, then type:
`Rscript bin/unscRamble.R path/to/Mapping_file.txt path/to/chr1_imputed.vcf.gz yes path/to/output_directory`

After the run, you will find all chromosomes in individual files (VCF format) in the desired folder.

# Questions/Comments/Feedback
For any question, comment, bug report or suggestion, please email `n.tesi@amsterdamumc.nl`.
