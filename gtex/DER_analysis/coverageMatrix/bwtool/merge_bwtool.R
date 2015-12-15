## Load required libraries
library('GenomicRanges')
library('BiocParallel')
library('devtools')

## Parallel environment to use
bp <- SnowParam(workers = 25, outfile = Sys.getenv('SGE_STDERR_PATH'))

## Load pheno data
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/pheno/pheno_missing_less_10.Rdata')

## Load bwtool results
tsv <- dir('/dcl01/leek/data/gtex_work/mean_cov', pattern = 'tsv', full.names = TRUE)
names(tsv) <- gsub('.mean.tsv', '', dir('/dcl01/leek/data/gtex_work/mean_cov/', pattern = 'tsv'))

## Find which samples have tsv files and match them that way
i <- match(as.character(pheno$Run), names(tsv))
if(any(is.na(i))) message(paste(Sys.time(), 'there are', sum(is.na(i)), 'tsv files missing'))

sampleNames <- gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', pheno$BigWigPath[!is.na(i)])

coverageMatrix <- bpmapply(function(tsvFile, sampleName) {
    message(paste(Sys.time(), 'reading file', tsvFile))
    res <- read.table(tsvFile, header = FALSE, colClasses = list(NULL, NULL, NULL, 'numeric'))
    colnames(res) <- sampleName
    return(as.matrix(res))
}, tsv[i[!is.na(i)]], sampleNames, BPPARAM = bp, SIMPLIFY = FALSE)
coverageMatrix <- do.call(cbind, coverageMatrix)

save(coverageMatrix, file = 'coverageMatrix.Rdata')


## Load regions
regs_raw <- read.table(tsv[1], header = FALSE, colClasses = list('character', 'numeric', 'numeric', NULL), col.names = c('chr', 'start', 'end', 'empty'))

chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
seqlen <- chrInfo$length[match(unique(regs_raw$chr), chrInfo$chr)]
names(seqlen) <- unique(regs_raw$chr)

regions <- GRanges(seqnames = regs_raw$chr, ranges = IRanges(regs_raw$start + 1, regs_raw$end), seqlengths = seqlen)
chr_n <- as.vector(table(regs_raw$chr))
names(chr_n) <- names(table(regs_raw$chr))
chr_n <- chr_n[match(unique(regs_raw$chr), names(chr_n))]
names(regions) <- paste(regs_raw$chr, unlist(sapply(chr_n, seq_len)), sep = '.')

save(regions, file = 'regions.Rdata')

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
