#!/usr/bin/Rscript

########################################################################
# Random set generation and scrambling                                 #
# The first operation consists in generating a set of random numbers   #
# These numbers will be assigned to each sample_ID.                    #
# The numbers will be shuffled in each chromosome.                     #
########################################################################

# Libraries
    library(argparse)
    library(data.table)
    library(stringr)
    args <- commandArgs(trailingOnly = FALSE)

# Functions
    # Function to take the VCF file (inputation input) and generate 22 PLINK files (1 per chromosome)
    function_split <- function(i, finp, out, const_FID){
            if (i == 1){ 
                if (const_FID == TRUE){
                    cmd = paste0("plink --vcf ", finp, " --real-ref-alleles --make-bed --out ", out, "/chrAll_input --threads 4 --const-fid"); system(cmd, ignore.stdout=T) 
                } else {
                    cmd = paste0("plink --vcf ", finp, " --real-ref-alleles --make-bed --out ", out, "/chrAll_input --threads 4"); system(cmd, ignore.stdout=T) 
                }
            }
            if (i == 23){ i = "X" }
            cmd = paste0("plink2 --bfile ", out, "/chrAll_input --chr ", i, " --real-ref-alleles --make-pgen --out ", out, "/chr", i, "_input --threads 4")
            system(cmd, ignore.stdout=T)
            cat(paste0("**** Done with chromosome ", i, "   \r"))
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
            if (sex_chromosomes == TRUE){
                for (chr in 1:23){ d[, paste0("CHR", chr)] <- sample(random_numbers) }
            } else {
                for (chr in 1:22){ d[, paste0("CHR", chr)] <- sample(random_numbers) }
            }

            # output the mapping file
            write.table(d, paste0(out, "/Mapping_file.txt"), quote=F, row.names=F, sep="\t")

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

            cat(paste0("**** Done with chromosome ", i, "    \r"))
    }

    # Function to merge back all files together
    function_merge <- function(out, sex_chromosomes){
            # first generate a list of files to merge to chr1...
            if (sex_chromosomes == TRUE){
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

# Main
# Parse arguments
    # Create parser
    parser <- ArgumentParser(description = "scRamble: a pipeline to scramble individual IDs in a VCF file across chromosomes, to make it virtually impossible to identify individuals in a dataset.")
    # Required arguments
        # Input VCF file
        parser$add_argument("--vcf", help = "Path to the input VCF file to be used for imputation.", type = "character", required = TRUE)
        # Output folder
        parser$add_argument("--out", help = "Path to the output folder where the results will be stored.", type = "character", required = TRUE)
        # Reference genome
        parser$add_argument("--ref", help="Reference genome build used.", type = "character", choices = c('hg19', 'hg38'), required = TRUE)
        # Sex chromosomes
        parser$add_argument("--sex", help="When present, sex chromosomes will also be scrambled.", action = "store_true", default = FALSE)
        # Flag for FID
        parser$add_argument("--constFID", help="When present, family IDs (FIDs) will be all set to default (Useful if sample name contains multiple '_').", action = "store_true", default = FALSE)

# Read arguments
    args = parser$parse_args()
    input_vcf = args$vcf
    output_folder = args$out
    reference_genome = args$ref
    sex_chromosomes = args$sex
    const_FID = args$constFID

# Print settings of the run
    cat("\n\n")
    cat("*** Welcome to scRamble ***\n")
    cat("** Your settings are:")
    cat(paste0("\n** Input VCF --> ", input_vcf))
    cat(paste0("\n** Output folder --> ", output_folder))
    cat(paste0("\n** Reference genome --> ", reference_genome))
    cat(paste0("\n** Consider sex chromosomes --> ", sex_chromosomes))
    cat(paste0("\n** Constant FID --> ", const_FID))
    cat("\n\n")

# Check if the input file exists
    if (!file.exists(input_vcf)){
        stop(paste0("** The input file ", input_vcf, " does not exist. Please check the path."))
    }

# Check if the output folder exists
    if (dir.exists(output_folder)){
        stop(paste0("** The output folder '", output_folder, "' already exists. We recommend using a fresh directory."), call. = FALSE)
    } else {
        success <- tryCatch({
            dir.create(output_folder, recursive = TRUE)
        }, error = function(e) {
            stop(paste0("** Failed to create output folder '", output_folder, "'. Error: ", conditionMessage(e)), call. = FALSE)
        })
        cat("** Created output folder:", output_folder, "\n")
    }

# First is the split of the main VCF into separate files per chromosome
    cat("** Start with splitting main VCF into chromosome-specific PLINK files.\n")
    if (sex_chromosomes == TRUE){
        res <- lapply(1:23, function_split, finp = input_vcf, out = output_folder, const_FID = const_FID)
    } else {
        res <- lapply(1:22, function_split, finp = input_vcf, out = output_folder, const_FID = const_FID)
    }
    cat("** Done with splitting main VCF into chromosome-specific PLINK files.\n\n")

# Second is the generation of the random numbers and their distribution across chromosomes
    cat("** Start with generation of random numbers.\n")
    key <- function_generate(output_folder, sex_chromosomes)
    cat("** Done with generation of random numbers.\n\n")

# Third we need to loop across chromosome and do the scrambling
    cat("** Start with scrambling IDs of each PLINK file.\n")
    if (sex_chromosomes == TRUE){ 
        res <- lapply(1:23, function_scramble, key = key, out = output_folder)
    } else {
        res <- lapply(1:22, function_scramble, key = key, out = output_folder)
    }
    cat("** Done with scrambling IDs of each PLINK file.\n\n")

# Fourth we need to re-merge all files together
    cat("** Start with merging.\n")
    res <- function_merge(output_folder, sex_chromosomes)
    cat("** Done with merging. A VCF file with scrambled individual IDs has been created.\n\n")

# Clean things up
    cat("** Removing redundant files.\n\n")
    try(system(paste0("rm ", output_folder, "/fileList_toMerge.txt")), silent=T)
    try(system(paste0("rm ", output_folder, "/*nosex")), silent=T)
    try(system(paste0("for chr in {1..22}; do rm ", output_folder, "/chr${chr}_*; done")), silent=T)
    try(system(paste0("rm ", output_folder, "/*log")))
    try(system(paste0("rm ", output_folder, "/*input.*")), silent=T)
    try(system(paste0("rm ", output_folder, "/*bed")), silent=T)
    try(system(paste0("rm ", output_folder, "/*bim")), silent=T)
    try(system(paste0("rm ", output_folder, "/*fam")), silent=T)

# Finally add 'chr' if Reference genome is GRCh38 and update the contig IDs
    if (toupper(reference_genome) == "HG38"){
        cat("** Adding chr in front of each chromosome as data is GRCh38.\n\n")
        cmd = paste0("zcat ", output_folder, "/chrAll_input_scrambled.vcf.gz | awk '{if($0 !~ /^#/) print \"chr\"$0; else print $0}' > ", output_folder, "/chrAll_input_scrambled.vcf")
        system(cmd)
        cmd = paste0("sed 's/contig=<ID=/contig=<ID=chr/g' ", output_folder, "/chrAll_input_scrambled.vcf | sed 's/ID=chr23/ID=chrX/g' > ", output_folder, "/chrAll_input_scrambled_updated.vcf")
        system(cmd)
        system(paste0("rm ", output_folder, "/chrAll_input_scrambled.*"))
        system(paste0("mv ", output_folder, "/chrAll_input_scrambled_updated.vcf ", output_folder, "/chrAll_input_scrambled.vcf"))
    } else {
        cmd = paste0('zcat ', output_folder, '/chrAll_input_scrambled.vcf.gz > ', output_folder, '/chrAll_input_scrambled.vcf')
        system(cmd)
    }

# Last thing is that if chrX is there it is called chr23 --> change the name now
    if (sex_chromosomes == TRUE){
        cat("** Fixing chromosome X as this is included in the data.\n\n")
        cmd = paste0("sed 's/chr23/chrX/g' ", output_folder, "/chrAll_input_scrambled.vcf > ", output_folder, "/chrAll_input_scrambled_updated.vcf")
        system(cmd)
        system(paste0("rm ", output_folder, "/chrAll_input_scrambled.*"))
        system(paste0("mv ", output_folder, "/chrAll_input_scrambled_updated.vcf ", output_folder, "/chrAll_input_scrambled.vcf"))
    } 

# Finally, sort, index and divide in chromosomes
    cat("** Sorting data.\n\n")
    cmd_sort = paste0("bcftools sort -Oz -o ", output_folder, "/chrAll_input_scrambled_sorted.vcf.gz ", output_folder, "/chrAll_input_scrambled.vcf > /dev/null 2>&1")
    system(cmd_sort)
    cat("** Indexing data.\n\n")
    cmd_index = paste0("bcftools index ", output_folder, "/chrAll_input_scrambled_sorted.vcf.gz ")
    system(cmd_index)
    cat("** Splitting chromosomes.\n\n")
    if (sex_chromosomes == TRUE){
        cmd_divide = paste0("for chr in {1..22} X; do bcftools view ", output_folder, "/chrAll_input_scrambled_sorted.vcf.gz --regions chr${chr} -Oz -o ", output_folder, "/chr${chr}_input_scrambled.vcf.gz; done")
    } else {
        cmd_divide = paste0("for chr in {1..22}; do bcftools view ", output_folder, "/chrAll_input_scrambled_sorted.vcf.gz --regions chr${chr} -Oz -o ", output_folder, "/chr${chr}_input_scrambled.vcf.gz; done")
    }
    system(cmd_divide)
    cat("** Cleaning data.\n\n")
    system(paste0("rm ", output_folder, "/chrAll*"))

# Final message
    cat("** Done. Goodbye!\n")

