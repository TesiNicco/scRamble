#!/usr/bin/Rscript

# Libraries
    library(argparse)
    library(data.table)
    library(stringr)
    args <- commandArgs(trailingOnly = FALSE)

# Functions
unscramble = function(i, key, filelist, out_path, single_chromosome){
    print(paste0("# Started with chromosome ", i))
    # Take chromosome of interest
    if (single_chromosome == TRUE){
        data_interest = filelist[1]
    } else {
        data_interest = filelist[i]
    }
    # grep index of FID old column
    index_fid = grep("FID", colnames(key))
    # grep index of the IID old column
    index_iid = grep("IID", colnames(key))
    # grep index of the FAID column
    index_faid = grep("FAID", colnames(key))
    # create dataframe with the old IDs (FID and IID) and the new IDs (FAID and chromosome-specific)
    key_interest = data.frame(key, check.names=F)[, c(index_fid, index_iid, index_faid, index_faid + i)]
    # rename columns
    colnames(key_interest) = c("real_family_name", "real_sample_name", "scrambled_name_family", "scrambled_name")
    # Extract samples names
    #header = unlist(strsplit(system(paste0("zcat ", data_interest, " | head -19 | tail -1"), inter=T), "\t"))
    header = unlist(strsplit(system(paste0("bcftools view -h ", data_interest, " | tail -1"), intern=T), "\t"))
    header = data.frame(scrambled = header[10:length(header)], order = seq(1, length(header[10:length(header)])))
    # Merge with mapping file
    header_merged = merge(header, key_interest, by.x = "scrambled", by.y = "scrambled_name")
    # Re-order -- just in case
    header_merged = header_merged[order(header_merged$order), ]
    # create new id column mergring family and sample name
    header_merged$real_name = paste(header_merged$real_family_name, header_merged$real_sample_name, sep="_")
    # Ok, now write the temporary mapping file
    write.table(header_merged$real_name, "tmp.txt", quote=F, row.names=F, col.names=F)
    # Now we have to use bcftools to update names
    cmd = paste0("bcftools reheader -s tmp.txt -o ", out_path, "/chr", i, ".dose.unscrambled.vcf.gz ", data_interest)  
    system(cmd)
    # Remove temporary mapping file
    system(paste0("rm tmp.txt"))
    return(paste("# Done with chromosome ", i))
}

# Main
# Parse arguments
    parser <- ArgumentParser(description='unscRamble: given scrambled VCF files and a mapping file, unscramble the VCF files to reconstruct the original individual genotypes.')
    # Required arguments
        # Mapping file
        parser$add_argument('--map', type='character', help='Mapping file with the scrambled and real IDs.', required=TRUE)
        # VCF file(s)
        parser$add_argument('--vcf', type='character', help='Input VCF file(s) with scrambled IDs. Imputation output is normally split by chromosome. You just need to indicate chr1 file (e.g chr1.imputed.vcf.gz)', required=TRUE)
        # Sex chromosomes
        parser$add_argument("--sex", help="When present, sex chromosomes will also be scrambled.", action = "store_true", default = FALSE)
        # Output directory
        parser$add_argument('--out', type='character', help='Output directory, where unscRambled VCF files will be placed.', required=TRUE)
        # Flag to indicate whether a single chromosome should be performed
        parser$add_argument("--single", help="When present, only a single chromosome will be processed.", action = "store_true", default = FALSE)

# Read arguments
    args = parser$parse_args()
    mapping_file = args$map
    input_vcf = args$vcf
    sex_chromosomes = args$sex
    output_folder = args$out
    single_chromosome = args$single

# Print settings of the run
    cat("\n\n")
    cat("*** Welcome to unscRamble ***\n")
    cat("** Your settings are:")
    cat(paste0("\n** Mapping file --> ", mapping_file))
    cat(paste0("\n** Input VCF --> ", input_vcf))
    cat(paste0("\n** Output folder --> ", output_folder))
    cat(paste0("\n** Consider sex chromosomes --> ", sex_chromosomes))
    cat(paste0("\n** Process a single chromosome --> ", single_chromosome))
    cat("\n\n")

# Check if the mapping file exists and in case read it
    if (!file.exists(mapping_file)){
        stop(paste0("** The mapping file ", mapping_file, " does not exist. Please check the path and try again."))
    }
    cat(paste0("** The mapping file exists.\n\n"))
    key <- fread(mapping_file, h=T)

# Check if the input VCF file exists and in case process it
    if (!file.exists(input_vcf)){
        stop(paste0("** The input VCF file ", input_vcf, " does not exist. Please check the path and try again."))
    }
    cat(paste0("** The input VCF file exists.\n\n"))
    if (single_chromosome == TRUE){
        cat("** Processing a single chromosome as requested.\n\n")
        filelist <- c(input_vcf)
        chrom_interest = str_replace_all(str_extract(input_vcf, "chr[0-9XY]+"), 'chr', '')
        chrom_interest = ifelse(chrom_interest == "X", 23, as.numeric(chrom_interest))
    } else {
        fnames <- str_split_fixed(input_vcf, "chr[0-22]", 2)
        if (sex_chromosomes == TRUE){
            filelist <- paste0(fnames[,1], "chr", c(seq(1,22), "X"), fnames[,2])
        } else {
            filelist <- paste0(fnames[,1], "chr", seq(1,22), fnames[,2])
        }
    }

# Check if output directory exists and in case create it
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

# Run for all chromosomes
    if (single_chromosome == TRUE){
        res <- unscramble(chrom_interest, key = key, filelist = filelist, out_path = output_folder, single_chromosome = TRUE)
    } else {
        if (sex_chromosomes == TRUE){
            res <- lapply(1:23, unscramble, key = key, filelist = filelist, out_path = output_folder, single_chromosome = FALSE)
        } else {
            res <- lapply(1:22, unscramble, key = key, filelist = filelist, out_path = output_folder, single_chromosome = FALSE)
        }
    }

