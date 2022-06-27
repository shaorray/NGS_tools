# R script for insertion size coverage

# Rui Shao, 2022, Jun

# --------------------------------------------- initialization ------------------------------------------------- #
args <- commandArgs(TRUE)

pkgs <- c("optparse", "foreach", "doParallel", "BiocManager")
suppressMessages(for (pkg in pkgs) {
  if (!require(pkg, character.only = T)) install.packages(pkg)
  require(pkg, quietly = T, character.only = T)
})

bioc_pkgs <- c("S4Vectors", "Rsamtools", "rtracklayer", "GenomeInfoDb", "IRanges",
               "GenomicRanges", "IRanges", "GenomicAlignments")
suppressMessages(for (pkg in bioc_pkgs) {
  if (!require(pkg, character.only = T)) BiocManager::install(pkg)
  require(pkg, quietly = T, character.only = T)
})

options("scipen" = 999, "digits" = 4)
options(warn = -1)

# --------------------------------------------- specifications ------------------------------------------------- #

option_list <- list(
  # bam input
  make_option(c("-b", "--bam"), action = "store", default = NULL, type = 'character',
              help = "Bam file, paired end."),
  
  make_option(c("-l", "--isize_lim"), action = "store", default = 2000, type = 'numeric',
              help = "Max reads insertion size [default %default]"),
  
  make_option(c("-q", "--min_mapq"), action = "store", default = 10, type = 'numeric',
              help = "Min reads map quality [default %default]"),
  
  make_option(c("-s", "--is_smooth"), action = "store", default = FALSE, type = 'logical',
              help = "Perform bin smooth? [default %default]"),
  
  make_option("--bin_length", action = "store", default = 1000, type = 'numeric',
              help = "Bin length, if is smooth [default %default]"),
  
  make_option(c("-t", "--thread"), action = "store", default = 1, type = 'numeric',
              help = "Number of CPU threads [default %default]"),
  # output
  make_option(c("-o", "--out"), action = "store", default = "frag_isize_coverage", type = 'character',
              help = "Output file name"),
  
  make_option("--quietly", action = "store", default = FALSE, type = 'logical',
              help = "Mute verbos and messages [default %default]")
  
)
arguments <- parse_args(OptionParser(option_list = option_list),
                        positional_arguments = TRUE)
opt <- arguments$options


# parse arguments
file_path <- opt$bam

isize_lim <- opt$isize_lim
min_mapq <- opt$min_mapq

is_smooth <- opt$is_smooth
bin_length <- opt$bin_length

file_name <- opt$out

is_queitly <- opt$quietly

if (!is_queitly) {
  message("\n----------------------------------------------\n")
  message(paste("Fragment isize coverage: ", date(), "\n"))
  message("----------------------------------------------\n")
  message(paste("BAM:\n\n", file_path, "\n"))
  message("----------------------------------------------\n")
  message(paste("Max insertion size:\n\n", isize_lim, "\n"))
  message("----------------------------------------------\n")
  message(paste("Min map quality:\n\n", min_mapq, "\n"))
  message("----------------------------------------------\n")
  message(paste("Perform bin smooth:\n\n", is_smooth, "\n"))
  message("----------------------------------------------\n")
  message(paste("Bin length:\n\n", bin_length, "\n"))
  message("----------------------------------------------\n")
  message(paste("Output file name:\n\n", file_name, "\n"))
  message("----------------------------------------------\n")
}

doParallel::registerDoParallel(cores = opt$thread)

# --------------------------------------------- functions ------------------------------------------------- #

bam_to_fragment_size_bw <- function(file_path,
                                    file_name,
                                    isize_lim = 2000,
                                    min_mapq = 10, 
                                    is_smooth = FALSE,
                                    bin_length = 1000) {
  
  sbp <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isFirstMateRead = TRUE,
                                                               isDuplicate = FALSE,
                                                               isSecondaryAlignment = FALSE),
                                 what = c("pos", "mpos", "isize", "mapq"),
                                 mapqFilter = min_mapq)
  gr <- GenomicAlignments::readGAlignments(file_path, param = sbp)
  
  gr <- gr[which(abs(rtracklayer::mcols(gr)$isize) < isize_lim)]
  
  gr_seqlengths <- seqlengths(gr)
  
  # convert fragment GRanges to coverage
  frag_gr <- GRanges(seqnames = seqnames(gr), 
                     IRanges(start = ifelse(strand(gr) == "+", 
                                            mcols(gr)$pos, mcols(gr)$mpos),
                             width = abs(mcols(gr)$isize)), 
                     strand = "*")
  frag_gr$isize <- width(frag_gr)
  seqlengths(frag_gr) <- gr_seqlengths
  frag_cov <- coverage(frag_gr, weight = as.numeric(frag_gr$isize)) / coverage(frag_gr) 
  
  # bin the coverage
  if (is_smooth) {
    chr_tile <- tileGenome(gr_seqlengths, tilewidth = bin_length, cut.last.tile.in.chrom = TRUE)
    chr_tile$score <- foreach (chr = seqlevels(chr_tile), .combine = c) %dopar% 
      {
        if (!chr %in% names(frag_cov)) return(rep(NA, sum(seqnames(chr_tile) == chr)))
        bin_mean(x = frag_cov[[chr]],
                 bin_width = bin_length,
                 stride = bin_length, 
                 add_name = F)
      }
    chr_tile <- chr_tile[!is.na(chr_tile$score)]
    chr_tile$score[which(chr_tile$score <= 0)] <- median(chr_tile$score[which(chr_tile$score > 0)])
    
    frag_cov <- coverage(chr_tile, weight = chr_tile$score)
  }
  
  # save coverage
  export.bw(frag_cov, paste0(file_name, '.bw')) # the empty regions are filled with NAs
}



bin_mean <- function(x, bin_width = NULL, stride = NULL, add_name = FALSE)
{
  # Args:
  #           x: a vector
  #   bin_width: the window for sum
  #      stride: bin step length
  #    add_name: bin names
  
  if (is.null(bin_width)) bin_width = 1000
  if (is.null(stride) | bin_width < stride) stride = bin_width
  
  x = as.numeric(x)
  len = length(x)
  
  # bin ranges
  starts = c(1, seq_len((len - 1) %/% stride) * stride )
  ends = starts + bin_width 
  ends[ends > len] = len
  
  # bin sums with stride
  x_is_na = is.na(x)
  x_0 = x
  x_0[x_is_na] = 0
  
  c_sums = cumsum(x_0)
  c_sums = c_sums[ends] - c_sums[starts] 
  
  effective_lens = cumsum(!x_is_na)
  effective_lens = effective_lens[ends] - effective_lens[starts] 
  
  c_means = c_sums / effective_lens
  
  if (add_name) names(c_means) = starts
  c_means
}

# --------------------------------------------------- run ------------------------------------------------------ #

bam_to_fragment_size_bw(file_path = file_path,
                        file_name = file_name,
                        isize_lim = isize_lim,
                        min_mapq = min_mapq, 
                        is_smooth = is_smooth,
                        bin_length = bin_length)
