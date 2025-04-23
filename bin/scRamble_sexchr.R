### Libraries
library(data.table)
library(stringr)
#library(ggplot2)
args = commandArgs(trailingOnly = TRUE)

### Functions
# Function to take the VCF file (inputation input) and generate 22 PLINK files (1 per chromosome)
function_split <- function(i, finp, out){
        if (i == 1){ cmd = paste0("plink --vcf ", finp, " --real-ref-alleles --make-bed --out ", out, "/chrAll_input --threads 4"); system(cmd, ignore.stdout=T) }
        if (i == 23){ i = "X" }
        cmd = paste0("plink2 --bfile ", out, "/chrAll_input --chr ", i, " --real-ref-alleles --make-pgen --out ", out, "/chr", i, "_input --threads 4")
        system(cmd, ignore.stdout=T)
        print(paste0("##        Done with chromosome ", i))
        return("Done with splitting")
}

# Function to generate and output the random numbers to use for scrambling
# This is important because these codes will be used to match things back
function_generate <- function(out, sex_chromosomes){
        # list all .psam files in the folder
        listFiles <- system(paste0("ls ", out, "/*psam"), intern=T)

        # open the first
        d <- fread(listFiles[1], h=T, check.names=F)

        # add family ID if missing
        if (!any(grepl("FID", colnames(d)))){
            d$FID = 0
        }

        # add new family ID as 0
        d$FAID = 0

        # generate as many random numbers as the number of samples
        set.seed(1234)
        random_numbers <- sample(x=seq(100000, 999999), size = nrow(d), replace=F)

        # for loop to assign a random number to each sample, and shuffle these in the chromosomes
        if (sex_chromosomes %in% c("yes", "YES")){
            for (chr in 1:23){ d[, paste0("CHR", chr)] <- sample(random_numbers) }
        } else {
            for (chr in 1:22){ d[, paste0("CHR", chr)] <- sample(random_numbers) }
        }

        # output the mapping file
        write.table(d, paste0(out, "/Mapping_file.txt"), quote=F, row.names=F, sep="\t")

        print("## Done with the generation of numbers")
        return(d)
}

# Function to scramble the individual IDs of ea chromosome
function_scramble <- function(i, key, out){
        # make sure key is a dataframe
        key = data.frame(key, check.names=F)

        # find the index of the chromosome
        index_chrom = which(colnames(key) == paste0('CHR', i))

        # find the index of the old FID column
        index_fid = grep('FID', colnames(key))

        # find the index of the IID column
        index_id = which(colnames(key) == 'IID')
        if (length(index_id) == 0){
            index_id = grep('#IID', colnames(key))
        }
        if (length(index_id) == 0){
            # stop the script with error
            stop("## Error: No IID or #IID column found in the mapping file. Please check the input file.")
        }

        # find the index of the new FID
        index_newfid = which(colnames(key) == 'FAID')

        # generate temporary file for the new list of keys
        write.table(key[, c(index_fid, index_id, index_newfid, index_chrom)], paste0(out, "/tmp_key.txt"), quote=F, row.names=F, sep="\t")

        if (i == 23){ i = "X" }

        # run /tudelft.net/staff-bulk/ewi/insy/DBL/svanderlee/programs/plink2 to update IDs
        cmd <- paste0("plink2 --pfile ", out, "/chr", i, "_input --update-ids ", out, "/tmp_key.txt --real-ref-alleles --make-bed --out ", out, "/chr", i, "_input_updated --threads 4")
        system(cmd, ignore.stdout=T)

        # remove key
        system(paste0("rm ", out, "/tmp_key.txt"))

        print(paste0("##        Done with chromosome ", i))
        return(paste0("Done with chromosome ", i))
}

# Function to merge back all files together
function_merge <- function(out, sex_chromosomes){
        # first generate a list of files to merge to chr1...
        if (sex_chromosomes %in% c("yes", "YES")){
            tmp <- paste0(out, "/chr", c(seq(2, 22), "X"), "_input_updated")
        } else {
            tmp <- paste0(out, "/chr", seq(2, 22), "_input_updated")
        }
        write.table(tmp, paste0(out, "/fileList_toMerge.txt"), quote=F, row.names=F, col.names=F)

        # second, the actual merge
        cmd = paste0("plink --bfile ", out, "/chr1_input_updated --merge-list ", out, "/fileList_toMerge.txt --real-ref-alleles --recode vcf-iid bgz --out ", out, "/chrAll_input_scrambled --threads 4")
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

# Print settings of the run
cat("\n\n")
cat("### Welcome to scRamble ###")
cat("\n")
cat("## Your settings are:")
cat(paste0("\n## INPUT VCF --> ", args[1]))
cat(paste0("\n## OUTPUT FOLDER --> ", args[2]))
cat(paste0("\n## CONSIDER SEX CHROMOSOMES --> ", args[3]))
cat(paste0("\n## REFERENCE GENOME --> ", args[4]))
cat("\n\n")

system(paste0("mkdir ", args[2]))

# First is the split of the main VCF into separate files per chromosome
print("## Start with splitting main VCF into chromosome-specific PLINK files.")
if (args[3] %in% c("yes", "YES")){ 
    res <- lapply(1:23, function_split, finp = args[1], out = args[2])
} else {
    res <- lapply(1:22, function_split, finp = args[1], out = args[2])
}
print("## Done with splitting main VCF into chromosome-specific PLINK files.")

# Second is the generation of the random numbers and their distribution across chromosomes
print("## Start with generation of random numbers.")
key <- function_generate(args[2], args[3])
print("## Done with generation of random numbers.")

# Third we need to loop across chromosome and do the scrambling
print("## Start with scrambling IDs of each PLINK file.")
if (args[3] %in% c("yes", "YES")){ 
    res <- lapply(1:23, function_scramble, key = key, out = args[2])
} else {
    res <- lapply(1:22, function_scramble, key = key, out = args[2])
}
print("## Done with scrambling IDs of each PLINK file.")

# Fourth we need to re-merge all files together
print("## Start with merging.")
res <- function_merge(args[2], args[3])
print("## Done with merging. A VCF file with scrambled individual IDs has been created.")

# Clean things up
print("## Removing redundant file.")
try(system(paste0("rm ", args[2], "/fileList_toMerge.txt")), silent=T)
try(system(paste0("rm ", args[2], "/*nosex")), silent=T)
try(system(paste0("for chr in {1..22}; do rm ", args[2], "/chr${chr}_*; done")), silent=T)
try(system(paste0("rm ", args[2], "/*log")))
try(system(paste0("rm ", args[2], "/*input.*")), silent=T)
try(system(paste0("rm ", args[2], "/*bed")), silent=T)
try(system(paste0("rm ", args[2], "/*bim")), silent=T)
try(system(paste0("rm ", args[2], "/*fam")), silent=T)

# Finally add 'chr' if Reference genome is GRCh38 and update the contig IDs
if (args[4] %in% c("GRCh38", "hg38", "HG38")){
    print("## Adding chr in front of each chromosome as data is GRCh38.")
    cmd = paste0("zcat ", args[2], "/chrAll_input_scrambled.vcf.gz | awk '{if($0 !~ /^#/) print \"chr\"$0; else print $0}' > ", args[2], "/chrAll_input_scrambled.vcf")
    system(cmd)
    cmd = paste0("sed 's/contig=<ID=/contig=<ID=chr/g' ", args[2], "/chrAll_input_scrambled.vcf | sed 's/ID=chr23/ID=chrX/g' > ", args[2], "/chrAll_input_scrambled_updated.vcf")
    system(cmd)
    system(paste0("rm ", args[2], "/chrAll_input_scrambled.*"))
    system(paste0("mv ", args[2], "/chrAll_input_scrambled_updated.vcf ", args[2], "/chrAll_input_scrambled.vcf"))
}

# Last thing is that if chrX is there it is called chr23 --> change the name now
if (args[3] %in% c("yes", "YES")){
    print("## Fixing chromosome X as this is included in the data.")
    cmd = paste0("sed 's/chr23/chrX/g' ", args[2], "/chrAll_input_scrambled.vcf > ", args[2], "/chrAll_input_scrambled_updated.vcf")
    system(cmd)
    system(paste0("rm ", args[2], "/chrAll_input_scrambled.*"))
    system(paste0("mv ", args[2], "/chrAll_input_scrambled_updated.vcf ", args[2], "/chrAll_input_scrambled.vcf"))
}

# Finally, sort, index and divide in chromosomes
print("## Sorting data.")
cmd_sort = paste0("bcftools sort -Oz -o ", args[2], "/chrAll_input_scrambled_sorted.vcf.gz ", args[2], "/chrAll_input_scrambled.vcf")
system(cmd_sort)
print("## Indexing data.")
cmd_index = paste0("bcftools index ", args[2], "/chrAll_input_scrambled_sorted.vcf.gz ")
system(cmd_index)
print("## Splitting chromosomes.")
cmd_divide = paste0("for chr in {1..22} X; do bcftools view ", args[2], "/chrAll_input_scrambled_sorted.vcf.gz --regions chr${chr} -Oz -o ", args[2], "/chr${chr}_input_scrambled.vcf.gz; done")
system(cmd_divide)
print("## Cleaning data.")
system(paste0("rm ", args[2], "/chrAll*"))
print("## Done. Goodbye!")

