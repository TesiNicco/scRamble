### Libraries
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

# Function to generate IBD and PCA plot of the actual correct files
function_PCA_IBD_preScramble <- function(out){
    # Generate a set of 20000 randomly selected common variants and extract these variants from main file
    system(paste0("mkdir ", out, "/1000Genome"), ignore.stdout=T)
    system(paste0("plink2 --pfile ", out, "/chrAll_input --maf 0.05 --make-pgen --out ", out, "/1000Genome/chrAll_input_commonOnly --threads 4"), ignore.stdout=T)
    system(paste0("grep -v '#' ", out, "/1000Genome/chrAll_input_commonOnly.pvar | sort -R | cut -f3,3 | uniq | head -20000 > ", out, "/1000Genome/random_20k_snps.txt"))

    # Do the same extraction in 1000Genome data -- then merge all 1000Genome data in 1 file
    system(paste0("for chr in {1..22}; do plink2 --bfile /tudelft.net/staff-bulk/ewi/insy/DBL/niccolo/databases/variants_DB/1kGenome/plink/chr${chr} --extract ", out, "/1000Genome/random_20k_snps.txt --make-bed --out ", out, "/1000Genome/chr${chr}_1000Genome --threads 4; done"), ignore.stdout=T)
    system(paste0("ls ", out, "/1000Genome/chr*1000Genome.bim | sed 's/.bim//g' | grep -v chr1_ > ", out, "/1000Genome/fileList_1000Genome.txt"))
    system(paste0("plink --bfile ", out, "/1000Genome/chr1_1000Genome --merge-list ", out, "/1000Genome/fileList_1000Genome.txt --make-bed --out ", out, "/1000Genome/chrAll_1000Genome_randomSNPs --threads 4"), ignore.stdout=T)

    # Then merge 1000Genome data with our actual data
    system(paste0("cut -f2,2 ", out, "/1000Genome/chrAll_1000Genome_randomSNPs.bim > ", out, "/1000Genome/snps_to_keep.txt"))
    system(paste0("plink2 --pfile ", out, "/chrAll_input --extract ", out, "/1000Genome/snps_to_keep.txt --make-bed --out ", out, "/1000Genome/chrAll_input_commonOnly_randomSNPs --threads 4"), ignore.stdout=T)
    system(paste0("rm ", out, "/1000Genome/chrAll_input_commonOnly.*"), ignore.stdout=T)
    system(paste0("plink --bfile ", out, "/1000Genome/chrAll_1000Genome_randomSNPs --bmerge ", out, "/1000Genome/chrAll_input_commonOnly_randomSNPs --make-bed --out ", out, "/1000Genome/chrAll_combined_randomSNPs --threads 4"), ignore.stdout=T)

    # Then do PCA
    print("##   Performing PCA")
    system(paste0("plink2 --bfile ", out, "/1000Genome/chrAll_combined_randomSNPs --pca --out ", out, "/1000Genome/combined_pca_preScramble --threads 4"), ignore.stdout=T)

    # Then do IBD -- this is without the 1000Genome data
    print("##   Calculating IBD")
    system(paste0("mkdir ", out, "/IBD"), ignore.stdout=T)
    system(paste0("plink --bfile ", out, "/1000Genome/chrAll_input_commonOnly_randomSNPs --genome --maf 0.2 --min 0.2 --out ", out, "/IBD/IBS_over0.2_relations --threads 4"), ignore.stdout=T)
}

# Function to generate IBD and PCA plot of the actual correct files
function_PCA_IBD_postScramble <- function(out){
    # Generate a set of 20000 randomly selected common variants and extract these variants from main file
    system(paste0("plink2 --bfile ", out, "/chrAll_input_scrambled --extract ", out, "/1000Genome/snps_to_keep.txt --threads 4 --make-bed --out ", out, "/1000Genome/chrAll_input_commonOnly_randomSNPs_scrambled"), ignore.stdout=T)

    # Merge with 1000Genome
    system(paste0("plink --bfile ", out, "/1000Genome/chrAll_1000Genome_randomSNPs --bmerge ", out, "/1000Genome/chrAll_input_commonOnly_randomSNPs_scrambled --make-bed --out ", out, "/1000Genome/chrAll_combined_randomSNPs_scrambled --threads 4"), ignore.stdout=T)

    # Then do PCA
    print("##   Performing PCA")
    system(paste0("plink2 --bfile ", out, "/1000Genome/chrAll_combined_randomSNPs_scrambled --pca --out ", out, "/1000Genome/combined_pca_postScramble --threads 4"), ignore.stdout=T)

    # Then do IBD -- this is without the 1000Genome data
    print("##   Calculating IBD")
    system(paste0("plink --bfile ", out, "/1000Genome/chrAll_input_commonOnly_randomSNPs_scrambled --genome --maf 0.2 --min 0.2 --out ", out, "/IBD/IBS_over0.2_relations_scrambled --threads 4"), ignore.stdout=T)
}

# Function to plot IBS
function_plotIBS <- function(finp){
        # Read input file
        ibd <- fread(finp, h=T)

        # Add colors
        ibd$col <- "darkolivegreen3"
        ibd$col[which(ibd$PI_HAT >0.3)] <- "grey40"
        ibd$col[which(ibd$PI_HAT >0.40)] <- "deepskyblue3"
        ibd$col[which(ibd$PI_HAT >0.80)] <- "red"
        
        # Sample points in the lower part of the plot
        if (nrow(ibd) > 0){
                rand_points <- data.frame(x = sample(1:nrow(ibd), 3000, replace=T), y = sample(seq(0, 0.2, 0.001), 3000, replace=T))
        } else {
                rand_points <- data.frame(x = sample(1:12000, 3000, replace=T), y = sample(seq(0, 0.2, 0.001), 3000, replace=T))
        }

        # Plot -- all pairs on the x -- PI_HAT on the y (PI_HAT measures overall IBD)
        plot(1:nrow(rand_points), rep(0, nrow(rand_points)), xlab="Pairs on individuals", ylab = "PI HAT (P(IBD=2) + 0.5*P(IBD=1))", xaxt="none", main="Identity by descent ~ pre-scramble",
                pch=16, col="white", ylim=c(0, 1.15), xlim=c(0, nrow(rand_points)))
        for (x in seq(0, 1.15, 0.115)){ abline(h=x, lwd=0.4, col="grey80") }
        for (x in seq(0, nrow(rand_points), nrow(rand_points)/10)){ abline(v=x, lwd=0.4, col="grey80") }
        if (nrow(ibd) >0){ points(1:nrow(ibd), ibd$PI_HAT, pch=16, col=alpha(ibd$col, 0.7)) }
        points(rand_points$x, rand_points$y, pch=16, col=alpha("darkolivegreen3", 0.7))
        abline(h=0.3, lty=2, col="red")
        legend("top", legend=c("Duplicate", "1 Degree relations", "More distant relations", "No relations"), pch=rep(16, 4), col=c("red", "deepskyblue3", "grey40", "darkolivegreen3"),
                bty="n", cex=1, pt.cex=1.50, ncol=2)
}

# Function to plot PCA
function_plotPCA <- function(finp){
        # Read input file
        pca <- fread(finp)

        # Read 1000Genome information
        k.fam <- read.table("~/bulkALL/niccolo/databases/variants_DB/1kGenome/popInfo/samples_description.txt", h=T)
      
        #merge phenotypes of 1K data with mds data
        fit.K.fam <- merge(pca, k.fam, by.x='IID', by.y='Sample', all.x=T)
        fit.K.fam$col <- ''
        fit.K.fam[fit.K.fam$Population == 'AFR',]$col <- 'navy'
        fit.K.fam[fit.K.fam$Population == 'EUR',]$col <- 'grey80'
        fit.K.fam[fit.K.fam$Population == 'SAS',]$col <- 'brown4'
        fit.K.fam[fit.K.fam$Population == 'EAS',]$col <- 'deepskyblue'
        fit.K.fam[fit.K.fam$Population == 'AMR',]$col <- 'gold1'
        fit.K.fam[fit.K.fam$Population == 'Han',]$col <- 'deepskyblue'
        fit.K.fam[fit.K.fam$Population == 'Chinese',]$col <- 'deepskyblue'
        fit.K.fam[fit.K.fam$Population == 'Vietnamese',]$col <- 'deepskyblue'
        fit.K.fam[fit.K.fam$Population == 'Mandinka',]$col <- 'navy'
        fit.K.fam[fit.K.fam$Population == 'SW',]$col <- 'navy'
        fit.K.fam[fit.K.fam$Population == 'Rican',]$col <- 'gold1'
        fit.K.fam[fit.K.fam$Population == 'Lankan',]$col <- 'brown4'
        fit.K.fam[fit.K.fam$col == '',]$col <- 'black'
        
        # plot
        plot(fit.K.fam$PC1[which(fit.K.fam$col != "black")], fit.K.fam$PC2[which(fit.K.fam$col != "black")], pch=16, col="white", cex=1.25, 
                xlab="PC1", ylab="PC2", xaxt="none", yaxt="none", main="Population stratification ~ pre-scramble")
        for (x in seq(min(fit.K.fam$PC1), max(fit.K.fam$PC1), (max(fit.K.fam$PC1)-min(fit.K.fam$PC1))/10)){ abline(v=x, lwd=0.4, col="grey80") }
        for (x in seq(min(fit.K.fam$PC2), max(fit.K.fam$PC2), (max(fit.K.fam$PC2)-min(fit.K.fam$PC2))/10)){ abline(h=x, lwd=0.4, col="grey80") }
        points(fit.K.fam$PC1[which(fit.K.fam$col != "black")], fit.K.fam$PC2[which(fit.K.fam$col != "black")], pch=16, col=alpha(fit.K.fam$col[which(fit.K.fam$col != "black")], 0.6), cex=1.25)
        points(fit.K.fam$PC1[which(fit.K.fam$col == "black")], fit.K.fam$PC2[which(fit.K.fam$col == "black")], pch=4, col=alpha(fit.K.fam$col[which(fit.K.fam$col == "black")], 0.6), cex=1)
        legend('topright', legend=c('Samples', 'EUR', 'AMR', 'AFR', 'SAS', 'EAS'), col=c('black', 'grey80', 'gold1', 'navy', 'brown4', 'deepskyblue2'), pch=c(4,16,16,16,16,16), ncol=2, bty="n")
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
        write.table(d, paste0(out, "/Mapping_file_input.txt"), quote=F, row.names=F, sep="\t")

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

# Function to check whether the scrambling has occurred
check_scramble <- function(out){
        print("## Checking whether scramble succedeed")

        # take randomly 1000 snps from the .bim file
        pvar <- fread(paste0(out, "/chrAll_input_scrambled.bim"), h=F)
        random_snps <- pvar[sample(1:nrow(pvar), 10000)]

        # extract these snps in 1 individual (1_33483 --> 838214 for chr1)
        sample_to_look <- data.frame(IID = c("838214", "1_33483"))
        write.table(random_snps$V2, paste0(out, "/tmp_check_scramble.txt"), quote=F, row.names=F, col.names=F)
        write.table(sample_to_look, paste0(out, "/tmp_check_scramble_sample.txt"), quote=F, row.names=F, col.names=F)
        cmd = paste0("plink2 --bfile ", out, "/chrAll_input_scrambled --extract ", out, "/tmp_check_scramble.txt --keep ", out, "/tmp_check_scramble_sample.txt --export A --out ", out, "/tmp_check_scramble_dosages.txt --threads 4")
        system(cmd, ignore.stdout=T)

        # for comparison do the same extract 
        cmd = (paste0("plink2 --pfile ", out, "/chrAll_input --extract ", out, "/tmp_check_scramble.txt --keep ", out, "/tmp_check_scramble_sample.txt --export A --out ", out, "/tmp_check_scramble_dosages_original.txt --threads 4"))
        system(cmd, ignore.stdout=T)

        # then read both dosages files
        dos_scramb <- fread(paste0(out, "/tmp_check_scramble_dosages.txt.raw"), h=T)
        dos_orig <- fread(paste0(out, "/tmp_check_scramble_dosages_original.txt.raw"), h=T)

        # parse them
        dos_scramb <- dos_scramb[, 7:ncol(dos_scramb)]
        dos_scramb_t <- as.data.frame(t(dos_scramb))
        dos_orig <- dos_orig[, 7:ncol(dos_orig)]
        dos_orig_t <-   as.data.frame(t(dos_orig))

        # match columns
        together <- merge(dos_scramb_t, dos_orig_t, by="row.names")
        colnames(together) <- c("SNP", "SC", "OR")
        together <- together[!is.na(together$SC),]
        together <- together[!is.na(together$OR),]

        # then do correlation
        print(paste0("Correlation between variants is ", cor(together$SC, together$OR)))

        # do the same on the subset of variants that are in chr1 -- cor should be 1 here
        chr_snps <- random_snps[which(random_snps$V1 == 1),]
        snp_info <- str_split_fixed(together$SNP, "_", 2)
        together$SNP_ID <- snp_info[, 1]
        chr1_match <- together[which(together$SNP_ID %in% chr_snps$V2),]
        print(paste0("Correlation between variants (this should be 1) is ", cor(chr1_match$SC, chr1_match$OR)))

        return("Check done!")
}

### Main
########################################################################
# Random set generation and scrambling
# The first operation consists in generating a set of random numbers
# These numbers will be assigned to each sample_ID.
# The numbers will be shuffled in each chromosome.
########################################################################

#system(paste0("mkdir ", args[2]))

# First is the split of the main VCF into separate files per chromosome
print("## Start with splitting main VCF into chromosome-specific PLINK files.")
#res <- lapply(1:22, function_split, finp = args[1], out = args[2]) 
print("## Done with splitting main VCF into chromosome-specific PLINK files.")

# Second is the generation of the random numbers and their distribution across chromosomes
print("## Start with generation of random numbers")
#key <- function_generate(args[2])
print("## Done with generation of random numbers")

# Third we need to loop across chromosome and do the scrambling
print("## Start with scrambling IDs of each PLINK file")
#res <- lapply(1:22, function_scramble, key = key, out = args[2])
print("## Done with scrambling IDs of each PLINK file")

# Fourth we need to re-merge all files together
print("## Start with merging.")
#res <- function_merge(args[2])
print("## Done with merging. A VCF file with scrambled individual IDs has been created")

# Maybe good to do some checks
# First, check whether the scrambling has occurred
#print(check_scramble(args[2]))

# Second, let's generate IBD and PCA plots at this stage with actual, correct data
# Only needed once -- takes a lot of time
print("## Generating PCA and IBD files for original file")
#res <- function_PCA_IBD_preScramble(args[2])
print("## Generating PCA and IBD files for scrambled file")
#res <- function_PCA_IBD_postScramble(args[2])

# Make PCA and IBD plots pre-scrambled
png(paste0(args[2], "/plots_preScramble.png"), height=8, width=15, res=300, units="in")
par(mfrow=c(1, 2))
function_plotIBS(paste0(args[2], "/IBD/IBS_over0.2_relations.genome"))
function_plotPCA(paste0(args[2], "/1000Genome/combined_pca_preScramble.eigenvec"))
dev.off()

# Make PCA and IBD plots post-scrambled
png(paste0(args[2], "/plots_postScramble.png"), height=8, width=15, res=300, units="in")
par(mfrow=c(1, 2))
function_plotIBS(paste0(args[2], "/IBD/IBS_over0.2_relations_scrambled.genome"))
function_plotPCA(paste0(args[2], "/1000Genome/combined_pca_postScramble.eigenvec"))
dev.off()

# Clean shit up
#try(system(paste0("rm ", args[2], "/fileList_toMerge.txt")), silent=T)
#try(system(paste0("rm ", args[2], "/*nosex")), silent=T)
#try(system(paste0("for chr in {1..22}; do rm ", args[2], "/chr${chr}_*; done")), silent=T)
#try(system(paste0("rm ", args[2], "/*input.*")), silent=T)
print("ok")
#try(system(paste0("rm ", args[2], "/*bed")), silent=T)
print("ok 2")
#try(system(paste0("rm ", args[2], "/*bim")), silent=T)
print("ok 3")
#try(system(paste0("rm ", args[2], "/*fam")), silent=T)
print("ok 4")
#try(system(paste0("rm ", args[2], "/1000Genome/*eigenval")), silent=T)
print("ok 5")
try(system(paste0("rm ", args[2], "/IBD/*nosex")), silent=T)
print("ok 6")
try(system(paste0("rm ", args[2], "/1000Genome/chr*combined*")), silent=T)
print("ok 7")
try(system(paste0("rm ", args[2], "/1000Genome/chr*input*")), silent=T)
print("top")
