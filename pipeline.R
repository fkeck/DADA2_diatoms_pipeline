library(dada2)

path <- "/media/francois/OS/dada_test/" # Your path to your sequences here

path_results <- file.path(path, "results")
if(!dir.exists(path_results)) dir.create(path_results)


#### REMOVAL OF PRIMERS ####
# This part is optional. If the primers have already been removed, skip it.
# This script skips all the verifications and checkup of the original tutorial
# If you are not sure about what you are doing, it is strongly advised to read
# the official DADA2 ITS Pipeline Workflow.
# https://benjjneb.github.io/dada2/ITS_workflow.html

fas_Fs_raw <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fas_Rs_raw <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

# This is our set of primers (from Vasselon et al. 2017 )
FWD <- c("AGGTGAAGTAAAAGGTTCWTACTTAAA",
         "AGGTGAAGTTAAAGGTTCWTAYTTAAA",
         "AGGTGAAACTAAAGGTTCWTACTTAAA")

REV <- c("CCTTCTAATTTACCWACWACTG",
         "CCTTCTAATTTACCWACAACAG")

path_cut <- file.path(path, "cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut)
fas_Fs_cut <- file.path(path_cut, basename(fas_Fs_raw))
fas_Rs_cut <- file.path(path_cut, basename(fas_Rs_raw))

R1_flags <- paste("-g", FWD, collapse = " ") 
R2_flags <- paste("-G", REV, collapse = " ") 

cutadapt <- "cutadapt" # Path to the executable
for(i in seq_along(fas_Fs_raw)) {
  cat("Processing", which(fas_Fs_raw == i), "/", length(fas_Fs_raw), "-------", i, "\n")
  system2(cutadapt, args = c(R1_flags, R2_flags,
                             "--discard-untrimmed",
                             "--max-n 0",
                             #paste0("-m ", 250-nchar(FWD)[1], ":", 250-nchar(REV)[1]), # Strong constraint: expected length
                             #paste0("-M ", 250-nchar(FWD)[1], ":", 250-nchar(REV)[1]), 
                             "-o", fas_Fs_cut[i], "-p", fas_Rs_cut[i],
                             fas_Fs_raw[i], fas_Rs_raw[i]))
}
out_1 <- cbind(ShortRead::qa(fas_Fs_raw)[["readCounts"]][,"read", drop = FALSE],
               ShortRead::qa(fas_Fs_cut)[["readCounts"]][,"read", drop = FALSE])

head(out_1)

###################################################

#### DADA2 Pipeline starts here ####

path_process <- path_cut # If you skipped primers removal, provide the path to your sequences here

fas_Fs_process <- sort(list.files(path_process, pattern = "_R1.fastq", full.names = TRUE))
fas_Rs_process <- sort(list.files(path_process, pattern = "_R2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fas_Fs_process), "_"), function(x) x[1])


#### Inspect read quality profiles ####
plotQualityProfile(fas_Fs_process[2])
plotQualityProfile(fas_Rs_process[2])

pdf(file.path(path_results, "Read_quality_profile_aggregated.pdf"))
  p <- plotQualityProfile(sample(fas_Fs_process, replace = FALSE,
                            size = ifelse(length(fas_Fs_process) < 100, length(fas_Fs_process), 100)),
                     aggregate = TRUE)
  p + ggplot2::labs(title = "Forward")
  p <- plotQualityProfile(sample(fas_Rs_process, replace = FALSE,
                            size = ifelse(length(fas_Rs_process) < 100, length(fas_Rs_process), 100)),
                     aggregate = TRUE)
  p + ggplot2::labs(title = "Reverse")
dev.off()


#### FILTER AND TRIM ####
fas_Fs_filtered <- file.path(path, "filtered", basename(fas_Fs_process))
fas_Rs_filtered <- file.path(path, "filtered", basename(fas_Rs_process))
all.equal(basename(fas_Fs_raw), basename(fas_Fs_filtered))

names(fas_Fs_filtered) <- sample.names
names(fas_Rs_filtered) <- sample.names

# Think twice before copying next command
# Check DADA2 official tutorial and the help of filterAndTrim function for details about arguments
out_2 <- filterAndTrim(fas_Fs_process, fas_Fs_filtered, fas_Rs_process, fas_Rs_filtered, truncLen = c(200, 170),
                     maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
head(out_2)


#### LEARN THE ERROR RATES ####
error_F <- learnErrors(fas_Fs_filtered, multithread = TRUE, randomize = TRUE)
error_R <- learnErrors(fas_Rs_filtered, multithread = TRUE, randomize = TRUE)

pdf(file.path(path_results, "Error_rates_learning.pdf"))
  p <- plotErrors(error_F, nominalQ = TRUE)
  p + ggplot2::labs(title = "Error Forward")
  p <- plotErrors(error_R, nominalQ = TRUE)
  p + ggplot2::labs(title = "Error Reverse")
dev.off()


#### DEREPLICATION, SAMPLE INFERENCE & MERGE PAIRED READS ####
merged_list <- vector("list", length(sample.names))
names(merged_list) <- sample.names

for(i in sample.names){
  cat("Processing", which(sample.names == i), "/", length(sample.names), "-------", i, "\n")
  derep_Fs <- derepFastq(fas_Fs_filtered[[i]], verbose = TRUE)
  derep_Rs <- derepFastq(fas_Rs_filtered[[i]], verbose = TRUE)
  dds_Fs <- dada(derep_Fs, err = error_F, multithread = TRUE, verbose = TRUE)
  dds_Rs <- dada(derep_Rs, err = error_R, multithread = TRUE, verbose = TRUE)
  merged_list[[i]] <- mergePairs(dds_Fs, derep_Fs, dds_Rs, derep_Rs, verbose = TRUE)
}

#### CONSTRUCT SEQUENCE TABLE ####
seqtab <- makeSequenceTable(merged_list)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# I like to save in CSV :)
# Note that if you need to reload these files to process them later with DADA2 functions
# (eg. merge several runs), it is better to save in .rds (see official tutorial).
write.csv(seqtab, file.path(path_results, "sequence_table.csv"))


#### REMOVE CHIMERA ####
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab_nochim)
write.csv(seqtab_nochim, file.path(path_results, "sequence_table_nochim.csv"))


#### TRACK READS THROUGH THE PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out_1, out_2, sapply(merged_list, getN), rowSums(seqtab_nochim))
colnames(track) <- c("raw", "noN", "cutadapt", "input", "filtered", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file.path(path_results, "track_reads.csv"))


#### ASSIGN TAXONOMY ####
# Here we download and use the Diat.barcode v7 pre-processed for DADA2.
# You can use your own local database if needed.
httr::GET("https://data.inra.fr/api/access/datafile/83898?gbrecs=true",
          httr::write_disk(tax_fas <- tempfile(fileext = ".gz")))
tax <- assignTaxonomy(seqtab_nochim, tax_fas, minBoot = 75,
                      taxLevels = c("Empire", "Kingdom", "Subkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                      outputBootstraps = TRUE, verbose = TRUE, multithread = TRUE)
write.csv(tax, file.path(path_results, "seq_nochim_tax.csv"))


httr::GET("https://data.inra.fr/api/access/datafile/83899?gbrecs=true",
          httr::write_disk(esp_fas <- tempfile(fileext = ".gz")))
exact_sp <- assignSpecies(seqtab_nochim, esp_fas)
write.csv(exact_sp, file.path(path_results, "seq_nochim_exact_sp.csv"))


