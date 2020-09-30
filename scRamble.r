library(data.table)
library(stringr)
library(ggplot2)
args = commandArgs(trailingOnly = TRUE)

### Functions
# Function to take the VCF file (inputation input) and generate 22 PLINK files (1 per chromosome)
function_split <- function(i, finp, out){
        if (i == 1){ cmd = paste0("plink2 --vcf ", finp, " --make-pgen --out ", out, "/chrAll_input --threads 4"); system(cmd, ignore.stdout=T) }
        cmd = paste0("plink2 --pfile ", out, "/chrAll_input --chr ", i, " --make-pgen --out ", out, "/chr", i, "_input --threads 4")
        system(cmd, ignore.stdout=T)

        print(paste0("##        Done with chromosome ", i))
        return("Done with splitting")
}

# Function to generate and output the random numbers to use for scrambling
# This is important because these codes will be used to match things back
function_generate <- function(out){
        # list all .psam files in the folder
        listFiles <- system(paste0("ls ", out, "/*psam"), intern=T)

        # open the first
        d <- fread(listFiles[1], h=T)

        # generate as many random numbers as the number of samples
        set.seed(1234)
        random_numbers <- sample(x=seq(100000, 999999), size = nrow(d), replace=F)

        # for loop to assign a random number to each sample, and shuffle these in the chromosomes
        for (chr in 1:22){ d[, paste0("CHR", chr)] <- sample(random_numbers) }

        # output the mapping file
        write.table(d, paste0(out, "/Mapping_file.txt"), quote=F, row.names=F, sep="\t")

        print("## Done with the generation of numbers")
        return(d)
}

# Function to scramble the individual IDs of ea chromosome
function_scramble <- function(i, key, out){
        # generate temporary file for the new list of keys
        write.table(key[, c(1, ..i+2)], paste0(out, "/tmp_key.txt"), quote=F, row.names=F, sep="\t")

        # run plink2 to update IDs
        cmd <- paste0("plink2 --pfile ", out, "/chr", i, "_input --update-ids ", out, "/tmp_key.txt --make-bed --out ", out, "/chr", i, "_input_updated --threads 4")
        system(cmd, ignore.stdout=T)

        # remove key
        system(paste0("rm ", out, "/tmp_key.txt"))

        print(paste0("##        Done with chromosome ", i))
        return(paste0("Done with chromosome ", i))
}

# Function to merge back all files together
function_merge <- function(out){
        # first generate a list of files to merge to chr1...
        tmp <- paste0(out, "/chr", seq(2, 22), "_input_updated")
        write.table(tmp, paste0(out, "/fileList_toMerge.txt"), quote=F, row.names=F, col.names=F)

        # second, the actual merge
        cmd = paste0("plink --bfile ", out, "/chr1_input_updated --merge-list ", out, "/fileList_toMerge.txt --recode vcf-iid bgz --out ", out, "/chrAll_input_scrambled --threads 4")
        system(cmd, ignore.stdout=T)

        return("Done merging")
}

### Main
########################################################################
# Random set generation and scrambling
# The first operation consists in generating a set of random numbers
# These numbers will be assigned to each sample_ID.
# The numbers will be shuffled in each chromosome.
########################################################################

system(paste0("mkdir ", args[2]))

# First is the split of the main VCF into separate files per chromosome
print("## Start with splitting main VCF into chromosome-specific PLINK files.")
res <- lapply(1:22, function_split, finp = args[1], out = args[2])
print("## Done with splitting main VCF into chromosome-specific PLINK files.")

# Second is the generation of the random numbers and their distribution across chromosomes
print("## Start with generation of random numbers")
key <- function_generate(args[2])
print("## Done with generation of random numbers")

# Third we need to loop across chromosome and do the scrambling
print("## Start with scrambling IDs of each PLINK file")
res <- lapply(1:22, function_scramble, key = key, out = args[2])
print("## Done with scrambling IDs of each PLINK file")

# Fourth we need to re-merge all files together
print("## Start with merging.")
res <- function_merge(args[2])
print("## Done with merging. A VCF file with scrambled individual IDs has been created")

# Clean shit up
try(system(paste0("rm ", args[2], "/fileList_toMerge.txt")), silent=T)
try(system(paste0("rm ", args[2], "/*nosex")), silent=T)
try(system(paste0("for chr in {1..22}; do rm ", args[2], "/chr${chr}_*; done")), silent=T)
try(system(paste0("rm ", args[2], "/*input.*")), silent=T)
try(system(paste0("rm ", args[2], "/*bed")), silent=T)
try(system(paste0("rm ", args[2], "/*bim")), silent=T)
try(system(paste0("rm ", args[2], "/*fam")), silent=T)