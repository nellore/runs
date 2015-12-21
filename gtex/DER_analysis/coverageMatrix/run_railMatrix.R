library('derfinder')
library('devtools')
library('BiocParallel')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'chr', 'c', 1, 'character', 'Chromosome in the following format: chr1, chrX, chrY',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}


## Options
cutoff <- 0.5

chrs <- opt$chr
## Get chr length
chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
chrlens <- chrInfo$length[chrInfo$chr %in% chrs]

load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/pheno/pheno_missing_less_10.Rdata')

summaryFiles <- '/dcl01/leek/data/gtex_work/gtex_mean_coverage.bw'
sampleFiles <- pheno$BigWigPath
names(sampleFiles) <- gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', sampleFiles)
sampleFiles <- head(sampleFiles, 500)

## Run railMatrix
regionMat <- railMatrix(chrs, summaryFiles, sampleFiles, L = 1, cutoff = cutoff, targetSize = 40e6 * 100, totalMapped = pheno$SumCoverage, file.cores = 1L, chunksize = 10000, verbose.load = FALSE, chrlens = chrlens)

## Save results
save(regionMat, file=paste0('regionMat-cut', cutoff, '-', opt$chr, '.Rdata'))

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
