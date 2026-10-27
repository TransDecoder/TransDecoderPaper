#!/usr/bin/env Rscript
# Simulate one chunk of the GM12878-like paired-end library with Polyester.
# usage: 04_polyester_chunk.R <chunk.fa> <chunk.counts.tsv> <outdir> <seed> <fraglen_mean> <fraglen_sd>
#   chunk.counts.tsv: transcript_id <tab> reads (rounded salmon NumReads), same order as chunk.fa
# Polyester settings: one sample, paired-end, 101 nt, normal fragment lengths (mean/sd from salmon's
# fld), error_model illumina5, strand_specific = TRUE (mate 1 = transcript sense; mates are swapped
# afterwards to give the RF/dUTP layout of the real library), fixed seed.
suppressMessages({library(polyester); library(Biostrings)})
# Two patches to Polyester 1.38.0 internals (applied in its namespace; nothing else changes):
# 1. generate_fragments: the normal fragment-length draw can give lengths <= 0, and subseq then fails
#    ("negative width"). Lengths below MIN_FRAG (50 nt) are redrawn (truncated normal).
# 2. add_platform_error: the illumina5 model covers read positions 1-100 only (plus a pos-0 identity
#    row), so 101-nt reads hit "NA in probability vector" at position 101. Position 101 uses the
#    position-100 error profile.
MIN_FRAG <- 50
patch <- function(fname, old, new) {
  f <- get(fname, envir = asNamespace("polyester"))
  src <- paste(deparse(f), collapse = "\n")
  stopifnot(lengths(regmatches(src, gregexpr(old, src, fixed = TRUE))) >= 1)
  f2 <- eval(parse(text = gsub(old, new, src, fixed = TRUE)))
  environment(f2) <- asNamespace("polyester")
  assignInNamespace(fname, f2, "polyester")
}
patch("generate_fragments", "fraglens = round(rnorm(L, mean = fraglen, sd = fragsd))",
      paste0("{ fraglens = round(rnorm(L, mean = fraglen, sd = fragsd)); ",
             "while (any(bad <- fraglens < ", MIN_FRAG, ")) ",
             "fraglens[bad] = round(rnorm(sum(bad), mean = fraglen, sd = fragsd)) }"))
patch("add_platform_error", "errMat = m[m$pos == pos, 2:6]", "errMat = m[m$pos == min(pos, max(m$pos)), 2:6]")
a <- commandArgs(trailingOnly = TRUE)
fa <- a[1]; ct <- a[2]; outdir <- a[3]; seed <- as.integer(a[4])
fl <- as.numeric(a[5]); fsd <- as.numeric(a[6])
tx <- readDNAStringSet(fa)
cnt <- read.delim(ct, header = FALSE, col.names = c("id", "n"), colClasses = c("character", "integer"))
stopifnot(identical(names(tx), cnt$id))
readmat <- matrix(cnt$n, ncol = 1)
cat(sprintf("chunk %s: %d transcripts, %d pairs, seed %d, fraglen %.2f sd %.2f\n",
            fa, length(tx), sum(cnt$n), seed, fl, fsd))
t0 <- Sys.time()
simulate_experiment_countmat(fasta = fa, readmat = readmat, outdir = outdir, paired = TRUE,
                             seed = seed, readlen = 101, fraglen = fl, fragsd = fsd,
                             distr = "normal", error_model = "illumina5", strand_specific = TRUE,
                             gzip = TRUE)
cat("elapsed", format(Sys.time() - t0), "\n")
cat("POLYESTER_CHUNK_DONE\n")
