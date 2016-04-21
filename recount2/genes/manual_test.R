## Usage:
# module load R/3.3
# mkdir -p logs
# Rscript manual_test.R > logs/manual_test_log.txt 2>&1

## Manual test for 1 sample
tsvFile <- '/dcl01/leek/data/recount2/coverage/SRR2064317.sum.tsv'
tsv <- read.table(tsvFile, header = FALSE, colClasses = list('character', 'numeric', 'numeric', 'numeric'), col.names = c('seqname', 'start', 'end', 'count'))

## Load bigwig
library('rtracklayer')
cov <- import('/dcl01/leek/data/sra/v2/batch_41/coverage_bigwigs/SRR2064317.bw', as = 'RleList')

## Load exons
load('ucsc-knowngene-hg38-exons.Rdata')

## Check the first gene
ex <- exons[[1]]

## Not the same as of https://github.com/nellore/runs/commit/e5f1d46517e7cc30ba9a49e9b7742066a43bf7c3
sum(Views(cov[[as.character(unique(seqnames(ex)))]], ranges(ex)))
sum(Views(cov[[as.character(unique(seqnames(ex)))]], start = tsv$start[1] + 1, end = tsv$end[1] + 1))
## bwtool looks at the start and end, not at the actual sub-ranges

## Fixed the previous issue at https://github.com/nellore/runs/commit/fa228e6e5627931286d37ee45d7b213f682cb828

## Group counts by gene
n <- elementNROWS(exons)
counts <- as.vector(tapply(tsv$count, rep(seq_len(length(n)), n), sum))

system.time(
    counts_bw <- sapply(exons, function(exon) {
        tryCatch(
            sum(sum(Views(cov[[as.character(unique(seqnames(exon)))]],
                ranges(exon)))),
            error = function(e) return(NA)
        )
    })
)

## Expect 0's for the NA cases
all(counts[is.na(counts_bw)] == 0)

## Expect same number for non NA-cases
all(counts[!is.na(counts_bw)] - counts_bw[!is.na(counts_bw)] == 0)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
