library('GenomicFeatures')
library('derfinder')
library('devtools')
library('TxDb.Hsapiens.UCSC.hg38.knownGene')

## Get data from Biomart
system.time(xx <- makeTxDbPackageFromBiomart(version = '0.99', maintainer = 'Leonardo Collado-Torres <lcollado@jhu.edu>', author = 'Leonardo Collado-Torres <lcollado@jhu.edu>', destDir = '~/'))

## Find sqlite file
sql_file <- dir(file.path('~', 'TxDb.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5', 'inst', 'extdata'), pattern = 'sqlite', full.names = TRUE)

## Load info
txdb <- loadDb(sql_file)

## Create genomic state for ENSEMBL
genomicState.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5 <- makeGenomicState(txdb=txdb, chrs= seqlevels(txdb), chrsStyle = NULL)
    
## Save for later use
save(genomicState.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5, 
    file='genomicState.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5.Rdata')
    
## Create genomic state for UCSC
genomicState.Hsapiens.UCSC.hg38.knownGene <- makeGenomicState(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, chrs = seqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene), chrsStyle = NULL)

## Save
save(genomicState.Hsapiens.UCSC.hg38.knownGene, file = 'genomicState.Hsapiens.UCSC.hg38.knownGene.Rdata')

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
