library('derfinder')
library('devtools')
library('BiocParallel')

## Options
cutoff <- 2.5

load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/pheno/pheno_missing_less_10.Rdata')

chrs <- paste0('chr', c(1:22, 'X', 'Y'))
summaryFiles <- rep('/dcl01/leek/data/gtex_work/gtex_mean_coverage.bw', 24)
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

## Check distribution
summary(counts$totalMapped)
summary(counts$totalMapped) / 1e6

## Run railMatrix
regionMat <- railMatrix(chrs, summaryFiles, sampleFiles, L = 1L, cutoff = cutoff, targetSize = 40e6, totalMapped = counts$totalMapped, BPPARAM.custom = MulticoreParam(workers = 24, outfile = Sys.getenv('SGE_STDERR_PATH')))

## Take into account that each sample had different lengths
for(chr in chrs) {
    regionMat[[chr]]$coverageMatrix <- regionMat[[chr]]$coverageMatrix / matrix(rep(pheno$avgLength, each = nrow(regionMat[[chr]]$coverageMatrix)), ncol = ncol(regionMat[[chr]]$coverageMatrix))
}

## Save results
save(regionMat, file=paste0('regionMat-cut', cutoff, '.Rdata'))

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
