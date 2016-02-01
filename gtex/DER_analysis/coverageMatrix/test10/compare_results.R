library('GenomicRanges')
library('devtools')

## Load R results
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/test10/regions-cut0.5.Rdata')
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/test10/coverageMatrix-cut0.5.Rdata')

## Load bwtool results
tsv <- dir('/dcl01/leek/data/gtex_work/mean_cov', pattern = 'tsv', full.names = TRUE)
names(tsv) <- gsub('.sum.tsv', '', dir('/dcl01/leek/data/gtex_work/mean_cov/', pattern = 'tsv'))
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/pheno/pheno_missing_less_10.Rdata')
i <- match(names(tsv), as.character(pheno$Run))

sampleNames <- gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', pheno$BigWigPath[i])

covMat <- mapply(function(tsvFile, sampleName) {
    message(paste(Sys.time(), 'reading file', tsvFile))
    res <- read.table(tsvFile, header = FALSE, colClasses = list(NULL, NULL, NULL, 'numeric'))
    colnames(res) <- sampleName
    return(as.matrix(res))
}, tsv, sampleNames, SIMPLIFY = FALSE)
covMat <- do.call(cbind, covMat)


regs_raw <- read.table(tsv[1], header = FALSE, colClasses = list('character', 'numeric', 'numeric', NULL), col.names = c('chr', 'start', 'end', 'empty'))

chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
seqlen <- chrInfo$length[match(unique(regs_raw$chr), chrInfo$chr)]
names(seqlen) <- unique(regs_raw$chr)

regs <- GRanges(seqnames = regs_raw$chr, ranges = IRanges(regs_raw$start + 1, regs_raw$end), seqlengths = seqlen)
chr_n <- as.vector(table(regs_raw$chr))
names(chr_n) <- names(table(regs_raw$chr))
chr_n <- chr_n[match(unique(regs_raw$chr), names(chr_n))]
names(regs) <- paste(regs_raw$chr, unlist(sapply(chr_n, seq_len)), sep = '.')


regs_ori <- regions
mcols(regs_ori) <- NULL

## Check that the regions are identical
identical(regs, regs_ori)


## Is the coverage matrix identical after sorting the variables to the same order?
covMat_sorted <- covMat[, colnames(coverageMatrix)]
identical(colnames(coverageMatrix), colnames(covMat_sorted))

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()

