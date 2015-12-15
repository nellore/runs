library('GenomicRanges')
library('devtools')

# Options
cutoff <- 0.5

chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
chrs <- chrInfo$chr

regionMat <- lapply(chrs, function(chr) {
    f <- paste0('regionMat-cut', cutoff, '-', chr, '.Rdata')
    if(!file.exists(f)) {
        message(paste(Sys.time(), 'failed to find', f))
        return(NULL)
    }
    load(f)
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

seqlen <- chrInfo$length[match(names(seqlengths(regions)), chrInfo$chr)]
names(seqlen) <- names(seqlengths(regions))
seqlengths(regions) <- seqlen

save(regions, file = paste0('regions-cut', cutoff, '.Rdata'))

## Extract coverage matrix
coverageMatrix <- lapply(regionMat, function(x) { x$coverageMatrix})
coverageMatrix <- do.call(rbind, coverageMatrix)
save(coverageMatrix, file = paste0('coverageMatrix-cut', cutoff, '.Rdata'))

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
