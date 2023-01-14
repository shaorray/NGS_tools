
# R script for ATAC-seq quality check with an aligned bam file

# Rui Shao, 2023, Jan

# ATAC QC features:
# - 1 total reads
# - 2 mapped percentage
# - 3 mapq > 30 percentage
# - 4 chrM
# - 5 FRiP (fraction of reads inside peak)
# - 6 TSS (+/- 2kb)

# ------------------------------------------------------------------------------
# initialization

args <- commandArgs(TRUE)

if (!require("optparse", character.only = TRUE)) install.packages("optparse")
suppressMessages(require("optparse", quietly = TRUE, character.only = TRUE))

options("scipen" = 999, "digits" = 4)
options(warn = -1)

# ------------------------------------------------------------------------------
# parse arguments

option_list <- list(
  # bam input
  make_option(c("-b", "--bam"), 
              action = "store", 
              default = NULL, 
              type = "character",
              help = "Bam file, paired end."),
  
  make_option(c("-p", "--peak"), 
              action = "store", 
              default = NULL, 
              type = "character",
              help = "ATAC called peaks [bed or narroPeak]."),
  
  make_option(c("-g", "--genome"), 
              action = "store",
              default = "mm10", 
              type = "character",
              help = "Reference genome in use [default %default]"),
  
  make_option(c("-l", "--lib_path"), 
              action = "store",
              default = NULL, 
              type = "character",
              help = "An optional R library path. [default %default]"),
  
  make_option(c("-d", "--file_dir"), 
              action = "store",
              default = ".", 
              type = "character",
              help = "Directory for ATAC report [default %default]")
  
)
arguments <- parse_args(OptionParser(option_list = option_list),
                        positional_arguments = TRUE)
opt <- arguments$options

# ------------------------------------------------------------------------------
# load required packages

# check and set lib path
lib_path <- opt$lib_path
if (!is.null(lib_path)) .libPaths(c(lib_path))

if (!require("pacman", character.only = TRUE)) install.packages("pacman")
pacman::p_load(IRanges,
               GenomeInfoDb,
               GenomicRanges,
               AnnotationDbi,
               GenomicFeatures, 
               Rsamtools, 
               rtracklayer)

# ------------------------------------------------------------------------------
# load options
bam_file <- opt$bam
peak_file <- opt$peak
file_dir <- opt$file_dir
which_genome <- opt$genome

message("\n----------------------------------------------\n")
message(paste(date(), "\n\nATAC QC", "\n"))
message("----------------------------------------------\n")
message(paste("Bam:\n", bam_file, "\n"))
message("----------------------------------------------\n")
message(paste("Peak:\n", peak_file, "\n"))
message("----------------------------------------------\n")
message(paste("Genome:\n", which_genome, "\n"))
message("----------------------------------------------\n")
message(paste("Out dir:\n", file_dir, "\n"))
message("----------------------------------------------\n")

# ------------------------------------------------------------------------------
# process arguments

# bam_file <- "/Volumes/dir/runs_Liu_lab/Spatial_ATAC_X101SC22082371-Z01-J154/bam/E9.5_10um-singlehole-Dia_DNA-1202-.bam"
if (!file.exists(paste0(bam_file, ".bai"))) indexBam(bam_file)
# bam_file <- BamFile(bam_file)


# peak_file <- "/Volumes/dir/runs_Liu_lab/Spatial_ATAC_X101SC22082371-Z01-J154/peak/E9.5_10um-singlehole-Dia_DNA-1202-_ATAC/E9.5_10um-singlehole-Dia_DNA-1202-_ATAC_peaks.narrowPeak"
peak.gr <- rtracklayer::import(peak_file)


# output directory of ATAC QC report
if (!dir.exists(file_dir)) dir.create(file_dir)
file_path <- paste0(file_dir, "/ATAC_report.txt")


# tss
if (which_genome == "mm10") {
  pacman::p_load(TxDb.Mmusculus.UCSC.mm10.ensGene)
  gene.gr <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm10.ensGene)
} else if (which_genome == "mm9") {
  pacman::p_load(TxDb.Mmusculus.UCSC.mm9.ensGene)
  gene.gr <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm9.ensGene)
} else if (which_genome == "hg19") {
  pacman::p_load(TxDb.Hsapiens.UCSC.hg19.knownGene)
  gene.gr <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
} else if (which_genome == "hg38") {
  pacman::p_load(TxDb.Hsapiens.UCSC.hg38.knownGene)
  gene.gr <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
} else {
  stop("Genome supported: mm9, mm10, hg19, hg38")
}

stopifnot(length(gene.gr) > 0)
gene.gr <- gene.gr[seqnames(gene.gr) %in% paste0("chr", c(1:100, "X", "Y", "M"))]

tss_2k.gr <- IRanges::trim(IRanges::promoters(sort(gene.gr), 
                                              upstream = 2000, 
                                              downstream = 2000))

tss_100.gr <- IRanges::trim(IRanges::promoters(sort(gene.gr), 
                                              upstream = 100, 
                                              downstream = 100))


# -------------------------------------------------------------------------------

if (TRUE) {
  
  if (!file.exists(file_path)) {
    write("bam_file peak_file which_genome date all_reads mapped_reads mapped_percent mapq30_percent chrM_percent  FRiP  FRiTSS_2k FRiTSS_100", 
          file = file_path)
  }
  
  out_list <- list()
  
  out_list$bam_file <- gsub(".*\\/(.*.bam)", "\\1", bam_file)
  out_list$peak_file <- gsub(".*\\/(.*)", "\\1", peak_file)
  out_list$which_genome <- which_genome
  out_list$date <- gsub(" ", "_", Sys.time())
  
  idx_stats <- Rsamtools::idxstatsBam(bam_file)
  
  # 1
  out_list$allR <- sum(idx_stats[, c("mapped", "unmapped")])
  
  # 2
  out_list$mappedR <- sum(idx_stats$mapped)
  out_list$mappedP <- 
    sum(idx_stats$mapped) / sum(idx_stats[, c("mapped", "unmapped")]) * 100
  
  # 3
  sbf <- Rsamtools::scanBamFlag(isPaired = TRUE,
                                isUnmappedQuery = FALSE,
                                hasUnmappedMate = FALSE,
                                isSecondaryAlignment = FALSE,
                                isDuplicate = FALSE,
                                isSupplementaryAlignment = FALSE)
  
  out_list$mapq <-
    Rsamtools::countBam(bam_file,
                        param = ScanBamParam(flag = sbf,
                                             mapqFilter = 30))[["records"]] /
    out_list$allR * 100
  
  # 4
  out_list$chrMP <- 
    idx_stats[idx_stats$seqnames == "chrM", "mapped"] / 
    sum(idx_stats$mapped) * 100
  
  # 5 FRiP
  out_list$FRiP <- 
    sum(countBam(bam_file,
                 param = ScanBamParam(flag = sbf,
                                      which = peak.gr))[["records"]]) /
    out_list$mappedR * 100
  
  # 6 TSS
  out_list$FRiTSS2k <-
    sum(countBam(bam_file,
                 param = ScanBamParam(flag = sbf,
                                      which = tss_2k.gr))[["records"]]) /
    out_list$mappedR * 100
  
  out_list$FRiTSS100 <-
    sum(countBam(bam_file,
                 param = ScanBamParam(flag = sbf,
                                      which = tss_2k.gr))[["records"]]) /
    out_list$mappedR * 100
  
  out <- unlist(out_list)[c("bam_file",  
                            "peak_file", 
                            "which_genome",
                            "date",  
                            "allR",
                            "mappedR", 
                            "mappedP",  
                            "mapq", 
                            "chrMP",
                            "FRiP",  
                            "FRiTSS2k",
                            "FRiTSS100")]
  out <- paste(out, collapse = "   ")
  write(out, file = file_path, append = TRUE)
}

message(paste("Processed file:", gsub(".*\\/(.*.bam)", "\\1", bam_file)))


