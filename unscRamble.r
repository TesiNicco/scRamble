# Libraries
library(data.table)
library(stringr)

# Functions
# Function to unscramble (this assumes that there is 1 PLINK2 file per chromosome)
function_unscramble <- function(i, key){
    # Read scrambled file (.psam)
    scramb <- fread("chr1_GSA12_input_scrambled.psam", h=T)

    # Generate temporary file for the new list of keys
    tmp <- key[, c(..i+2, 1)]
    tmp <- tmp[order(tmp$CHR1),]
    write.table(tmp, paste0("chr", i, "_tmp_mapping.txt"), quote=F, row.names=F, col.names=F, sep="\t")

    # run plink2 to update IDs
    cmd <- paste0("plink2 --pfile chr", i, "_GSA12_input_scrambled --update-ids chr", i, "_tmp_mapping.txt --make-pgen --out chr", i, "_GSA12_input_updated")
    system(cmd)

    # remove key
    system(paste0("rm chr", i, "_tmp_mapping.txt"))

    return(paste0("Done with chromosome ", i))
}

# Function to check whether everything went fine
function_check <- function(){
    print("## Checking whether UNscramble succedeed")

    # take randomly 1000 snps from the .bim file (chr1 for easyness)
    pvar <- fread("chr1_GSA12_input_updated.pvar", h=F)
    random_snps <- pvar[sample(1:nrow(pvar), 10000)]

    # extract these snps in 1 individual (1_71457 --> 215409 for chr1)
    sample_to_look <- data.frame(IID = "1_71457")
    write.table(random_snps$V3, "tmp_check_unscramble.txt", quote=F, row.names=F, col.names=F)
    write.table(sample_to_look, "tmp_check_unscramble_sample.txt", quote=F, row.names=F, col.names=F)
    cmd = "plink2 --pfile chr1_GSA12_input_updated --extract tmp_check_unscramble.txt --keep tmp_check_unscramble_sample.txt --export A --out tmp_check_scramble_dosages.txt"
    system(cmd)

    # for comparison do the same extract 
    cmd = paste0("plink2 --vcf ", finp, " --chr 1 --keep tmp_check_unscramble_sample.txt --make-pgen --out tmp_check_chr1_original")
    system(cmd)
    cmd = paste0("plink2 --pfile tmp_check_chr1_original --extract tmp_check_unscramble.txt --export A --out tmp_check_scramble_dosages_original.txt")
    system(cmd)

    # read files back in
    unsc <- fread("tmp_check_scramble_dosages.txt.raw", h=T)
    orig <- fread("tmp_check_scramble_dosages_original.txt.raw", h=T)

    # parse them
    unsc <- unsc[, 7:ncol(unsc)]
    dos_unsc_t <- as.data.frame(t(unsc))
    orig <- orig[, 7:ncol(orig)]
    dos_orig_t <-   as.data.frame(t(orig))

    # match columns
    together <- merge(dos_unsc_t, dos_orig_t, by="row.names")
    colnames(together) <- c("SNP", "UNSC", "OR")
    together <- together[!is.na(together$UNSC),]
    together <- together[!is.na(together$OR),]

    # then do correlation
    print(paste0("Correlation between variants is ", cor(together$UNSC, together$OR)))

    return("## Done with check.")
}

# Main
# Read mapping file
key <- fread("Mapping_file_GSA12.txt", h=T)

# Run function to unscramble
res <- lapply(1:22, function_unscramble, key = key)

# Check whether everything went fine
print(function_check())

# Clean up shit
system("rm tmp*")