## Load required libraries
library('derfinder')
library('GenomicRanges')
library('devtools')
library('rtracklayer')

chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
cutoff <- 0.5

## Find the regions for all the chromosomes given a specific cutoff
getRegs <- function(chr, chrlen, cutoff) {
    message(paste(Sys.time(), 'loading meanCov for', chr))
    meanCov <- loadCoverage('/dcl01/leek/data/gtex_work/gtex_mean_coverage.bw', chr, chrlen = chrlen)
    
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

regions <- mapply(getRegs, chrInfo$chr, chrInfo$length, MoreArgs = list(cutoff = cutoff))

regions <- unlist(GRangesList(regions))

seqlen <- chrInfo$length[match(names(seqlengths(regions)), chrInfo$chr)]
names(seqlen) <- names(seqlengths(regions))
seqlengths(regions) <- seqlen

## Save regions for the given cutoff
save(regions, file = paste0('regions-cut', cutoff, '.Rdata'))
export(regions, con = paste0('regions-cut', cutoff, '.bed'), format='BED')

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
session_info()



