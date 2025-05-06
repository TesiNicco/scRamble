## scRamble
**scRamble** is an R pipeline designed to scramble chromosomes across individuals to protect genetic data privacy. This is crucial before submitting data to non-European imputation panels due to EU privacy laws.

## In a nutshell
The GDPR regulations prohibit sharing genetic data outside the EU, posing challenges for collaboration and data processing on non-EU servers. Genotype imputation, often hosted in the US (https://imputation.biodatacatalyst.nhlbi.nih.gov/), requires uploading genetic data, conflicting with GDPR. scRamble addresses this by shuffling chromosomes across individuals, making genome reconstruction impossible without a key. It generates random sample identifiers and shuffles data while retaining a local key for reassembly. Post-shuffling, the data can be used for genotype imputation without compromising quality. After imputation, scRamble reorders chromosomes to their original state.

## Installation
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

## How to use
By running:  
```console
./bin/scRamble.R -h
./bin/unscRamble.R -h
```
will display the help message. `scRamble` parameters are:
- `--vcf` (**Mandatory**): path to genotype data in VCF format. This is often the input data for genotype imputation.
- `--out` (**Mandatory**): path to the output folder where results will be stored. Please use a new folder (`scRamble` will create the folder if not present).
- `--ref` (**Mandatory**): reference genome build used. Choices are hg19 or hg38.
- `--sex` (*Optional, flag*): flag to indicate whether sex chromosomes should be considered or not.

`unscRamble` parameters are;
- `--map` (**Mandatory**): path to the Mapping file, created with `scRamble`.
- `--vcf` (**Mandatory**): path to the imputed VCF files. Normally, these come split by chromosome. You just need to indicate the path to chr1 file (e.g, chr1.imputed.vcf.gz).
- `--out` (**Mandatory**): path to the output folder where results will be stored. Please use a new folder (`unscRamble` will create the folder if not present).
- `--sex` (*Optional, flag*): flag to indicate whether sex chromosomes should be considered or not.

## Examples
Run `scRamble`:
```sh
scRamble.R --vcf example_data/chrAll_1000Genome_sample10K_hg19.vcf.gz --out example_data/example_hg19 --ref hg19 [--sex (add the flag to consider X chromosome)]
scRamble.R --vcf example_data/chrAll_1000Genome_sample10K_hg38.vcf.gz --out example_data/example_hg38 --ref hg38 [--sex (add the flag to consider X chromosome)]
```
Depending on the number of samples and number of variants, this can take several minutes.  

Run `unscRamble`:
```sh
unscRamble.R --map example_data/example_hg19/Mapping_file.txt --vcf example_data/example_hg19/chr1_input_scrambled.vcf.gz --out example_data/example_hg19/unscRamble [--sex (add the flag to consider X chromosome)]
unscRamble.R --map example_data/example_hg38/Mapping_file.txt --vcf example_data/example_hg38/chr1_input_scrambled.vcf.gz --out example_data/example_hg38/unscRamble [--sex (add the flag to consider X chromosome)]
```
Depending on the number of samples and number of variants, this can take several minutes.  

## Questions/Comments/Feedback
For questions, comments, bug reports, or suggestions, email `n.tesi@amsterdamumc.nl`.