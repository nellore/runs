## Load required libraries
library('GenomicRanges')
library('BiocParallel')
library('Hmisc')
library('devtools')

## Parallel environment to use
bp <- SnowParam(workers = 20, outfile = Sys.getenv('SGE_STDERR_PATH'))

## Load pheno data
load('/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata')
metadata <- metadata[!is.na(metadata$bigwig_path), ]

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

## Split regions by chr, then into chunks
chr_n <- table(seqnames(regions))
chunks <- unlist(sapply(chr_n, function(n) {
    if(n <= 1e4) {
        return(rep('_chunk1', n))
    } else {
        x <- cut2(seq_len(n), m = 1e4)
        levels(x) <- paste0('_chunk', seq_len(length(levels(x))))
        return(as.character(x))
    }
}))
stopifnot(length(chunks) == length(regions))
chunk_grp <- paste0(seqnames(regions), chunks)

regions_split <- split(regions, chunk_grp)

## Load data for each chunk and save it

mapply(function(regions_to_subset, chunk_name) {
    message(paste(Sys.time(), 'processing', chunk_name))
    
    ## Identify which regions to keep
    regions_keep <- which(countOverlaps(regions, regions_to_subset) > 0)
    stopifnot(length(regions_keep) == length(regions_to_subset))
    
    ## Save actual subset of regions used
    regions_subset <- regions[regions_keep]
    save(regions_subset, file = paste0('regions_', chunk_name, '-cut0.5.Rdata'))
    
    ## Compute chunk coverage matrix
    coverageMatrix <- bpmapply(function(tsvFile, sampleName, regions_keep) {
        message(paste(Sys.time(), 'reading file', tsvFile))
        res <- tryCatch(
            read.table(tsvFile, header = FALSE, colClasses = list(NULL, NULL,
                NULL, 'numeric'))[regions_keep, , drop = FALSE],
                error = function(e) {
                    data.frame(empty = rep(0, length(regions_keep)))
                }
        )
        colnames(res) <- sampleName
        res <- as.matrix(res)
        return(res)
        }, tsv[i[!is.na(i)]], sampleNames, BPPARAM = bp, SIMPLIFY = FALSE,
            MoreArgs = list(regions_keep = regions_keep)
    )
    names(coverageMatrix) <- NULL
    coverageMatrix <- do.call(cbind, coverageMatrix)

    save(coverageMatrix, file = paste0('coverageMatrix-cut0.5-', chunk_name,
        '.Rdata'))
    
    ## Finish
    return(invisible(NULL))
})

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
