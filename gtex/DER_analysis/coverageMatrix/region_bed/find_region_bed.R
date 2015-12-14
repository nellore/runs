## Load required libraries
library('derfinder')
library('GenomicRanges')
library('BiocParallel')
library('devtools')
library('rtracklayer')

chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
cuts <- 0.5

## Parallel environment to use
bp <- SerialParam()

## Find the regions for all the chromosomes given a specific cutoff
getRegChr <- function(cutoff, chr, meanCov) {
    suppressPackageStartupMessages(library('derfinder'))
    suppressPackageStartupMessages(library('IRanges'))
    message(paste(Sys.time(), 'processing', chr, 'with cutoff', cutoff))
    regs <- findRegions(position = Rle(TRUE, length(meanCov$coverage[[1]])), fstats = meanCov$coverage[[1]], chr = chr, maxClusterGap = 300L, cutoff = cutoff, verbose = FALSE)
    if(is.null(regs)) {
        message(paste(Sys.time(), 'found no regions for', chr, 'with cutoff', cutoff))
        regs <- GRanges()
    }
    else if(!is(regs, 'GRanges')) {
        message(paste(Sys.time(), 'processing failed for', chr, 'with cutoff', cutoff))
        regs <- GRanges()
    }
    names(regs) <- NULL
    return(regs)
}
getRegs <- function(chr, chrlen, cutoffs, param) {
    message(paste(Sys.time(), 'loading meanCov for', chr))
    meanCov <- loadCoverage('/dcl01/leek/data/gtex_work/gtex_mean_coverage.bw', chr, chrlen = chrlen)
    res <- bplapply(cutoffs, getRegChr, chr = chr, meanCov = meanCov, BPPARAM = param)
    names(res) <- as.character(cutoffs)
    return(res)
}

region_cuts_raw <- mapply(getRegs, chrInfo$chr, chrInfo$length, MoreArgs = list(cutoffs = cuts, param = bp))
names(region_cuts_raw) <- chrInfo$chr

region_cuts <- lapply(as.character(cuts), function(cutoff, chromosomes = chrInfo$chr) {
    regs <- lapply(chromosomes, function(chr) { region_cuts_raw[[chr]][[cutoff]] })
    names(regs) <- chromosomes
    regs <- unlist(GRangesList(regs))
    return(regs)
})
names(region_cuts) <- as.character(cuts)

## Save regions for the given cutoff
regions <- region_cuts[[1]]
save(regions, file = paste0('regions-cut', cuts, '.Rdata'))
export(regions, con = paste0('regions-cut', cuts, '.bed'), format='BED')

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
session_info()



