## Usage
# module load R/3.3
# mkdir -p logs
# Rscript bplength.R > logs/bplength_log.txt 2>&1

## Create genes GRanges object with total (reduced) exon bp length
## Also save how the exons are related, for speeding up the tsv -> count matrix step
library('GenomicRanges')
library('S4Vectors')

## Load data
load('ucsc-knowngene-hg38-genes.Rdata')
load('ucsc-knowngene-hg38-exons.Rdata')

## Add length of reduced exons by gene
genes$bp_length <- sum(width(exons))
save(genes, file = 'ucsc-knowngene-hg38-genes-bp-length.Rdata')

## See length difference summary
summary(width(genes) - genes$bp_length)

## Group counts by gene
n <- elementNROWS(exons)
count_groups <- rep(seq_len(length(n)), n)
save(count_groups, file = 'count_groups.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
