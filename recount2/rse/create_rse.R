## Load libraries
library('getopt')
suppressPackageStartupMessages(library('GenomicRanges'))
suppressPackageStartupMessages(library('SummarizedExperiment'))

## Specify parameters
spec <- matrix(c(
    'project', 'p', 1, 'character', 'Project ID',
	'metadata', 'm', 1, 'character', 'Metadata file name',
	'projectid', 'i', 1, 'character', 'Project ID',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) {
    opt <- list(project = 'sra', 'metadata' = 'metadata_sra.Rdata',
        projectid = 'DRP000499')
    ## Largest one, to find memory needed
    opt <- list(project = 'sra', 'metadata' = 'metadata_sra.Rdata',
        projectid = 'SRP025982')
}

## Create output dir
dir.create(paste0('rse_', opt$project), showWarnings = FALSE)

## Load GRanges and metadata
load('/dcl01/leek/data/gtex_work/runs/recount2/genes/ucsc-knowngene-hg38-genes-bp-length.Rdata')
load('/dcl01/leek/data/gtex_work/runs/recount2/genes/ucsc-knowngene-hg38-exons.Rdata')
load('/dcl01/leek/data/gtex_work/runs/recount2/genes/count_groups.Rdata')
load(opt$metadata)

## Subset to project of interest
metadata <- subset(metadata, project == opt$projectid)
if(nrow(metadata) == 0) stop(paste('Invalid project id', opt$projectid))

## Subset to only use the samples that have tsv files
metadata <- metadata[!is.na(metadata$tsv_path), ]
if(nrow(metadata) == 0) stop(paste('No samples have bwtool tsv files for project', opt$projectid))
rownames(metadata) <- NULL

## Create output dir for the project
outdir <- paste0('rse_', opt$project, '/', opt$projectid)
dir.create(outdir, showWarnings = FALSE)

## Read counts from bwtool tsv output files
counts <- mapply(function(tsvFile, sampleName) {
    message(paste(Sys.time(), 'reading file', tsvFile))
    res <- read.table(tsvFile, header = FALSE, colClasses = list(NULL, NULL, NULL, 'numeric'))
    colnames(res) <- sampleName
    return(as.matrix(res))
}, metadata$tsv_path, metadata$run, SIMPLIFY = FALSE)
counts <- do.call(cbind, counts)

## Memory used by counts
print('Memory used by exon counts')
print(object.size(counts), units = 'Mb')

## Save exon counts
message(paste(Sys.time(), 'writing file', file.path(outdir, 'counts_exon.tsv')))
write.table(counts, file = file.path(outdir, 'counts_exon.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)
system(paste('gzip', file.path(outdir, 'counts_exon.tsv')))

## Remove bigwig and tsv file paths
metadata_clean <- metadata[, !colnames(metadata) %in% c('bigwig_path',
    'tsv_path')]

## Create exon level rse
exons_all <- unlist(exons)
rse_exon <- SummarizedExperiment(assays = list('counts' = counts),
    colData = DataFrame(metadata_clean), rowRanges = exons_all)
message(paste(Sys.time(), 'writing file', file.path(outdir, 'rse_exon.Rdata')))
save(rse_exon, file = file.path(outdir, 'rse_exon.Rdata'))

## Summarize counts at gene level
counts_gene <- lapply(split(as.data.frame(counts), count_groups), colSums)
counts_gene <- do.call(rbind, counts_gene)
rownames(counts_gene) <- names(genes)

## Memory used by counts at gene level
print('Memory used by gene counts')
print(object.size(counts_gene), units = 'Mb')

## Save gene counts
message(paste(Sys.time(), 'writing file', file.path(outdir, 'counts_gene.tsv')))
write.table(counts_gene, file = file.path(outdir, 'counts_gene.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)
system(paste('gzip', file.path(outdir, 'counts_gene.tsv')))

## Create gene level rse
rse_gene <- SummarizedExperiment(assays = list('counts' = counts_gene),
    colData = DataFrame(metadata_clean), rowRanges = genes)
message(paste(Sys.time(), 'writing file', file.path(outdir, 'rse_gene.Rdata')))
save(rse_gene, file = file.path(outdir, 'rse_gene.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
