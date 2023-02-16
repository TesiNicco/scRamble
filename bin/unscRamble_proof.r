# Libraries
library(data.table)
library(stringr)
args = commandArgs(trailingOnly = TRUE)

# Functions
# Function to unscramble (this assumes that there is 1 PLINK2 file per chromosome)
function_unscramble <- function(i, key, out_dir, conv_filelist){
    # Grep file related to chromosome i
    fint <- conv_filelist[grep(paste0("chr", i, "_"), conv_filelist)]
    
    # Read scrambled file (.psam)
    scramb <- fread(paste0(fint, ".psam"), h=T)

    # Generate temporary file for the new list of keys
    tmp <- key[, c(..i+2, 1)]
    tmp <- tmp[order(tmp[,1]),]
    write.table(tmp, paste0(out_dir, "/tmp_mapping.txt"), quote=F, row.names=F, col.names=F, sep="\t")

    # Generate output name
    fout = paste0(fint, "_unscrambled")

    # run plink2 to update IDs
    cmd <- paste0("plink2 --pfile ", fint, " --update-ids ", out_dir, "/tmp_mapping.txt --make-pgen --out ", fout)
    system(cmd)

    # remove key
    system(paste0("rm ", out_dir, "/tmp_mapping.txt"))

    # remove scrambled file
    system(paste0("rm ", fint, ".*"))

    return(paste0("Done with chromosome ", i))
}

# Function to convert VCF files into plink2 files
function_convert <- function(i, filelist, out_dir){
    # Get file to work on
    fn <- filelist[i]

    # Change output name
    fname_tmp = unlist(strsplit(fn, "/"))
    fname_tmp2 = str_split_fixed(fname_tmp[length(fname_tmp)], ".vcf.gz", 2)
    
    # The convert
    cmd = paste0("plink2 --vcf ", fn, " --make-pgen --out ", out_dir, "/", fname_tmp2[,1])
    system(cmd)

    return("Done")
}

# Function to check whether everything went fine
function_check <- function(){
    print("## Checking whether UNscramble succedeed")

    # take randomly 1000 snps from the .bim file (chr16 for easyness)
    pvar <- fread("out_unscambled/chr16_imputed_unscrambled.pvar", h=F)
    psam <- fread("out_unscambled/chr16_imputed_unscrambled.psam", h=F)
    random_snps <- pvar[sample(1:nrow(pvar), 1000)]
    random_samples <- psam[sample(1:nrow(psam), 10)]
    write.table(random_snps$V3, "out_unscambled/tmp_check_unscramble.txt", quote=F, row.names=F, col.names=F)
    write.table(random_samples$V1, "out_unscambled/tmp_check_unscramble_sample.txt", quote=F, row.names=F, col.names=F)
    
    # extract info from unscrambled and original sample
    cmd = "plink2 --pfile out_unscambled/chr16_imputed_unscrambled --extract out_unscambled/tmp_check_unscramble.txt --keep out_unscambled/tmp_check_unscramble_sample.txt --export A --out out_unscambled/tmp_unscrambled_genotypes"
    system(cmd)
    cmd = "plink2 --pfile data_unscRamble/for_proof/chrAll_original --extract out_unscambled/tmp_check_unscramble.txt --keep out_unscambled/tmp_check_unscramble_sample.txt --export A --out out_unscambled/tmp_original_genotypes"
    system(cmd)

    # read files back in
    unsc <- fread("out_unscambled/tmp_unscrambled_genotypes.raw", h=T)
    orig <- fread("out_unscambled/tmp_original_genotypes.raw", h=T)

    # sort by name
    unsc <- unsc[order(unsc$IID),]
    orig <- orig[order(orig$IID),]

    # parse them
    unsc <- unsc[, 7:ncol(unsc)]
    dos_unsc_t <- as.data.frame(t(unsc))
    orig <- orig[, 7:ncol(orig)]
    dos_orig_t <-   as.data.frame(t(orig))

    # match columns
    together <- merge(dos_unsc_t, dos_orig_t, by="row.names")
    colnames(together) <- c("SNPS", paste0("UNSC_", seq(1, 10)), paste0("ORIG_", seq(1, 10)))

    # loop and for each sample do correlation
    correlations <- c()
    for (i in 2:11){
        cs <- cor(na.exclude(together[, i]), na.exclude(together[, i+10]))
        correlations <- c(correlations, cs)
    }

    return("## Done with check.")
}

# Main
# Read mapping file
key <- fread(args[1], h=T)

# Identify files name
fnames <- str_split_fixed(args[2], "chr[0-22]", 2)
filelist <- paste0(fnames[,1], "chr", seq(1,22), fnames[,2])

# Create output directory first
system(paste0("mkdir ", args[3]))

# Then we need to convert to plink2 as most of the times this is VCF file
res <- lapply(1:22, function_convert, filelist = filelist, out_dir = args[3])

# Get all converted files
conv_filelist <- system(paste0("ls ", args[3], "/*log | sed 's/.log//g'"), intern=T)

# Run function to unscramble
res <- lapply(1:22, function_unscramble, key = key, out_dir = args[3], conv_filelist)

# At this point, let's take a sample of 1000 variants in 100 individuals and see correlation
print(function_check)