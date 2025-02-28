# New way to unscramble without converting to plink -- do everything from the VCF file

# Libraries
library(data.table)
library(stringr)
args = commandArgs(trailingOnly = TRUE)

# Functions
unscramble = function(i, key, filelist, out_path){
    print(paste0("# Started with chromosome ", i))
    # Take chromosome of interest
    data_interest = filelist[i]
    key_interest = data.frame(key)[, c(1, i+2)]
    colnames(key_interest) = c("real_name", "scrambled_name")
    # Extract samples names
    #header = unlist(strsplit(system(paste0("zcat ", data_interest, " | head -19 | tail -1"), inter=T), "\t"))
    header = unlist(strsplit(system(paste0("bcftools view -h ", data_interest, " | tail -1"), intern=T), "\t"))
    header = data.frame(scrambled = header[10:length(header)], order = seq(1, length(header[10:length(header)])))
    # Merge with mapping file
    header_merged = merge(header, key_interest, by.x = "scrambled", by.y = "scrambled_name")
    # Re-order -- just in case
    header_merged = header_merged[order(header_merged$order), ]
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
# Read mapping file
key <- fread(args[1], h=T)

# Identify files name
fnames <- str_split_fixed(args[2], "chr[0-22]", 2)
if (args[3] %in% c("yes", "YES")){
    filelist <- paste0(fnames[,1], "chr", c(seq(1,22), "X"), fnames[,2])
} else {
    filelist <- paste0(fnames[,1], "chr", seq(1,22), fnames[,2])
}

# Create output directory first
system(paste0("mkdir ", args[4]))
out_path = args[4]

# Run for all chromosomes
if (args[3] %in% c("yes", "YES")){
    res <- lapply(1:23, unscramble, key = key, filelist = filelist, out_path = out_path)
} else {
    res <- lapply(1:22, unscramble, key = key, filelist = filelist, out_path = out_path)
}

