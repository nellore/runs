## Extract exons for each gene
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library('rtracklayer')

## Get genes
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
save(genes, file = 'ucsc-knowngene-hg38-genes.Rdata')

## Get Exons
exons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')
save(exons, file ='ucsc-knowngene-hg38-exons.Rdata')

## Export exons as a BED file
export(exons, con = 'ucsc-knowngene-hg38.bed', format='BED')


## Print some info
length(genes)
length(exons)
length(exons) - length(genes)

## Names in exons not in genes
genes_miss <- names(exons)[!names(exons) %in% names(genes)]
length(genes_miss)
save(genes_miss, file = 'genes_miss.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
