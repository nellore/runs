## Extract exons for each gene
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library('rtracklayer')

## Keep only canonical chrs + chrM
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- keepSeqlevels(txdb, paste0('chr', c(1:22, 'X', 'Y', 'M')))

## Get genes
genes <- genes(txdb)
save(genes, file = 'ucsc-knowngene-hg38-genes.Rdata')

## Get Exons
exons <- exonsBy(txdb, by = 'gene')

## Reduce by gene exons won't be overlapping
exons <- reduce(exons)

save(exons, file ='ucsc-knowngene-hg38-exons.Rdata')

## Identify weird cases
strand <- strand(exons)
strand_run <- sapply(strand, nrun)
seq <- seqnames(exons)
seq_run <- sapply(seq, nrun)

## Weird cases
exons[which(strand_run > 1 & seq_run > 1)]
exons_exclude <- exons[strand_run > 1 | seq_run > 1]

## Export exons as a BED file
export(exons, con = 'ucsc-knowngene-hg38.bed', format='BED')

## Print some info
length(genes)
length(exons)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
