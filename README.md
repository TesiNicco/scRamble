# scRamble
**scRamble** is an R pipeline designed to scramble chromosomes across individuals to protect genetic data privacy. This is crucial before submitting data to non-European imputation panels due to EU privacy laws.

# In a nutshell
The GDPR regulations prohibit sharing genetic data outside the EU, posing challenges for collaboration and data processing on non-EU servers. Genotype imputation, often hosted in the US (https://imputation.biodatacatalyst.nhlbi.nih.gov/), requires uploading genetic data, conflicting with GDPR. scRamble addresses this by shuffling chromosomes across individuals, making genome reconstruction impossible without a key. It generates random sample identifiers and shuffles data while retaining a local key for reassembly. Post-shuffling, the data can be used for genotype imputation without compromising quality. After imputation, scRamble reorders chromosomes to their original state.

# Installation
Clone the `scRamble` package:
```sh
git clone https://github.com/TesiNicco/scRamble.git
```
Requirements:
1. `R` with `data.table` and `stringr` libraries (`install.packages('data.table')`, `install.packages('stringr')`)
2. `bcftools` (https://samtools.github.io/bcftools/bcftools.html)
3. `plink` (version 1.9 and 2, https://www.cog-genomics.org/plink/)

Install packages with Conda using the provided `.yml` file:
```sh
conda env create --name scRamble --file=scRamble/bin/scRamble_env.yml
conda activate scRamble
```
You may need to install `R` packages manually.

# Input data
Provide a single VCF file with all SNPs and individuals, and specify the output directory for scrambled genotypes.

# How to run
Run scRamble:
```sh
Rscript path/to/scRamble/bin/scRamble.R [path/to/input.vcf.gz] [path/to/output_directory_name] [consider sex chromosome (yes/no)] [reference genome (hg19/hg38)]
```
Example:
```sh
Rscript bin/scRamble_sexchr.R chrAll.input.vcf.gz scrambled_genotypes yes hg38
```
Ensure data alignment with the reference genome. The process may take several minutes.

# Output scrambled data
The output folder will contain:
1. Scrambled VCF files, one per chromosome, sorted and indexed.
2. A mapping file for unscrambling (`Mapping_file.txt`).

# Reconstruct original data after imputation
To unscramble and reconstruct original genotypes post-imputation:
```sh
Rscript path/to/bin/unscRamble.R [path/to/Mapping_file.txt] [path/to/chr1_imputed.vcf.gz] [consider sex chromosome (yes/no)] [path/to/output_directory]
```
Example:
```sh
Rscript bin/unscRamble.R path/to/Mapping_file.txt path/to/chr1_imputed.vcf.gz yes path/to/output_directory
```
The output folder will contain individual VCF files for each chromosome.

# Questions/Comments/Feedback
For questions, comments, bug reports, or suggestions, email `n.tesi@amsterdamumc.nl`.