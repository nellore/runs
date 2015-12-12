library('GenomicRanges')
library('devtools')

# Options
cutoff <- 3

chrs <- paste0('chr', c(1:22, 'X', 'Y'))

regionMat <- lapply(chrs, function(chr) {
    load(paste0('regionMat-cut', cutoff, '-', chr, '.Rdata'))
    res <- regionMat
    return(res)
})

## Merge
regionMat <- do.call(c, regionMat)

## Save merged results
save(regionMat, file=paste0('regionMat-cut', cutoff, '.Rdata'))

## Extract regions
regions <- lapply(regionMat, function(x) { x$regions})
regions <- unlist(GRangesList(regions))
save(regions, file = 'regions.Rdata')

## Extract coverage matrix
coverageMatrix <- lapply(regionMat, function(x) { x$coverageMatrix})
coverageMatrix <- do.call(rbind, coverageMatrix)
save(coverageMatrix, file = 'coverageMatrix.Rdata')

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
