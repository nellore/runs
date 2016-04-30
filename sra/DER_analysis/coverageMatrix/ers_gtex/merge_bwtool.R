## Load required libraries
library('GenomicRanges')
library('BiocParallel')
library('devtools')

## Parallel environment to use
bp <- MulticoreParam(workers = 20, outfile = Sys.getenv('SGE_STDERR_PATH'))

## Load pheno data
load('/dcl01/leek/data/gtex_work/runs/recount2/metadata/metadata_sra.Rdata')

## Load bwtool results
tsv <- dir('/dcl01/leek/data/sra_work/mean_cov_ers_gtex', pattern = 'tsv', full.names = TRUE)
names(tsv) <- gsub('.sum.tsv', '', dir('/dcl01/leek/data/sra_work/mean_cov_ers_gtex', pattern = 'tsv'))

## Find which samples have tsv files and match them that way
i <- match(metadata$run, names(tsv))
if(any(is.na(i))) message(paste(Sys.time(), 'there are', sum(is.na(i)), 'tsv files missing'))
    
## Load regions from the first tsv file
regs_raw <- read.table(tsv[1], header = FALSE, colClasses = list('character', 'numeric', 'numeric', NULL), col.names = c('chr', 'start', 'end', 'empty'))

## Get chr length info
chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
seqlen <- chrInfo$length[match(unique(regs_raw$chr), chrInfo$chr)]
names(seqlen) <- unique(regs_raw$chr)

## Define regions GRanges object
regions <- GRanges(seqnames = regs_raw$chr, ranges = IRanges(regs_raw$start + 1, regs_raw$end), seqlengths = seqlen)
chr_n <- as.vector(table(regs_raw$chr))
names(chr_n) <- names(table(regs_raw$chr))
chr_n <- chr_n[match(unique(regs_raw$chr), names(chr_n))]
names(regions) <- paste(regs_raw$chr, unlist(sapply(chr_n, seq_len)), sep = '.')

save(regions, file = 'regions-cut0.5.Rdata')

## Define sample names
sampleNames <- metadata$run[!is.na(i)]

## Load regions to subset
#load('some_file.Rdata') ## Assume it has a GRanges object called regions_to_subset

## Identify which regions to keep
#regions_keep <- which(countOverlaps(regions, regions_to_subset) > 0)
regions_keep <- seq_len(length(regions))

## Save actual subset of regions used
regions_subset <- regions[regions_keep]
save(regions_subset, file = 'regions_subset-cut0.5.Rdata')

## Compute overall coverage matrix, but subsetted to the regions of interest
coverageMatrix <- bpmapply(function(tsvFile, sampleName) {
    message(paste(Sys.time(), 'reading file', tsvFile))
    res <- tryCatch(
        read.table(tsvFile, header = FALSE, colClasses = list(NULL, NULL, NULL,
            'numeric'))[regions_keep, ],
            error = function(e) {
                data.frame(empty = rep(0, length(regions_keep)))
            }
    )
    colnames(res) <- sampleName
    res <- DataFrame(res) ## To compress memory by sample
    return(res)
}, tsv[i[!is.na(i)]], sampleNames, BPPARAM = bp, SIMPLIFY = FALSE)
coverageMatrix <- do.call(cbind, coverageMatrix)

save(coverageMatrix, file = 'coverageMatrix-cut0.5.Rdata')

## Rename the big matrix since we want all the small ones (by chr) to be called 'coverageMatrix'
covMat <- coverageMatrix

## Export coverage matrices by chr
s <- split(seq_len(length(regions)), seqnames(regions))
for(chr in names(s)) {
    message(paste(Sys.time(), 'creating coverage matrix for', chr))
    coverageMatrix <- covMat[s[[chr]], ]
    size <- as.vector(object.size(coverageMatrix)) / (1024^2)
    if(size < 1024) {
        message(paste(Sys.time(), 'coverage matrix for', chr, 'uses', round(size, 2), 'Mb in RAM'))
    } else {
        message(paste(Sys.time(), 'coverage matrix for', chr, 'uses', round(size / 1024, 2), 'Gb in RAM'))
    }
    save(coverageMatrix, file = paste0('coverageMatrix-cut0.5-', chr, '.Rdata'))
}

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
