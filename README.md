############################################################
############################################################
This is scRambleeeee

A pipeline written in R to scramble chromosomes across a set of individuals to avoid concerns about privacy of genetic data.
This should be done prior to submittion to non-European imputation panels due to privacy laws in the European union.

To run scRamble, move into the scRamble directory, then type:
Rscript bin/scRamble_sexchr.R [input.vcf.gz] [output_directory_name] [consider sex chromosome (yes/no)] [reference genome (hg19/hg38)]

For example, assuming that (i) the input file is called /path/to/chrAll.input.vcf.gz, (ii) you want the output files to 
be placed in folder path/to/scrambled_genotypes, you are interested in sex chromosomes and your genome version is hg38 then type:
Rscript bin/scRamble_sexchr.R path/to/chrAll.input.vcf.gz path/to/scrambled_genotypes yes hg38

After the run, you will find in the desired folder three output files:
1. the scrambled file in VCF format [chrAll_input_scrambled.vcf.gz]
2. mapping file to be used for unscrambling [Mapping_file.txt]

As it is virtually impossible to reconstruct the entire individual genotype without the key, these file can be shared and used for imputation.

###########################################################
###########################################################

After imputation is performed, you want to unscramble the outputs and reconstruct the original genotypes.

To do so, move into the scRamble directory and type:
Rscript bin/unscRamble.R [mapping_file.txt] [chr1_imputed.vcf.gz] [output_directory]

For example, assuming that (i) mapping file is called Mapping_file.txt, (ii) the imputed file, usually already divided in chromosomes, 
is called chr1_imputed.vcf.gz (make sure to specificy chromosome 1 here) and (iii) the output directory is out_unscrambled, then type:
Rscript bin/unscRamble.R path/to/Mapping_file.txt path/to/chr1_imputed.vcf.gz path/to/output_directory

After the run, you will find all chromosomes in individual files (PLINK2 format) in the desired folder.

###########################################################
###########################################################

For any question, comment, bug report or suggestion, please email n.tesi@amsterdamumc.nl

