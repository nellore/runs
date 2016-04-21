## Extract exons for each gene
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library('rtracklayer')

## Get genes with default option single.strand.genes.only = TRUE
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
save(genes, file = 'ucsc-knowngene-hg38-genes.Rdata')

## Get Exons
exons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')

## Keep only exons for gene ids that we selected previously
exons <- exons[names(exons) %in% names(genes)]

## Reduce exons by gene so the exons won't be overlapping each other inside a gene
exons <- reduce(exons)
save(exons, file ='ucsc-knowngene-hg38-exons.Rdata')

## Export exons as a BED file
export(unlist(exons), con = 'ucsc-knowngene-hg38.bed', format='BED')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
