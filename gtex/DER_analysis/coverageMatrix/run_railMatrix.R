library('derfinder')
library('devtools')
library('BiocParallel')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'chrnum', 'c', 1, 'character', 'Chromosome in the following format: 1, X, Y',
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
cutoff <- 2.5

load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/pheno/pheno_missing_less_10.Rdata')

chrs <- paste0('chr', opt$chrnum)
summaryFiles <- '/dcl01/leek/data/gtex_work/gtex_mean_coverage.bw'
sampleFiles <- pheno$BigWigPath

## Find count files
counts_files <- file.path(dir('/dcl01/leek/data/gtex', pattern = 'batch', full.names = TRUE), 'cross_sample_results', 'counts.tsv.gz')
names(counts_files) <- dir('/dcl01/leek/data/gtex', pattern = 'batch')

## Read in counts info
counts <- lapply(counts_files, function(i) {
    read.table(i, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
})
counts <- do.call(rbind, counts)
counts$totalMapped <- as.integer(sapply(strsplit(counts$total.mapped.reads, ','), '[[', 1))

## Match files to counts
map <- match(gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', sampleFiles), counts$X)
counts <- counts[map, ]

## Run railMatrix
regionMat <- railMatrix(chrs, summaryFiles, sampleFiles, L = pheno$avgLength, cutoff = cutoff, targetSize = 40e6, totalMapped = counts$totalMapped, file.cores = 10)

## Save results
save(regionMat, file=paste0('regionMat-cut', cutoff, '-chr', opt$chrnum, '.Rdata'))

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
