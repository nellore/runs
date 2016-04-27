## Prepare metadata
# module load R/3.3
# mkdir -p logs
# Rscript create_metadata.R -p "sra" > logs/create_metadata_sra_log.txt 2>&1
# Rscript create_metadata.R -p "gtex" > logs/create_metadata_gtex_log.txt 2>&1

library('getopt')

## Specify parameters
spec <- matrix(c(
	'project', 'p', 1, 'character', 'Project ID: sra or gtex',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

if(opt$project == 'sra') {
    ## Load sra info
    metadata <- read.table('/dcl01/leek/data/gtex_work/runs/sra/v2/recount2_metadata.tsv',
    header = TRUE, sep = '\t', stringsAsFactors = FALSE
    )
    colnames(metadata) <- tolower(colnames(metadata))
    colnames(metadata) <- gsub('\\.', '_', colnames(metadata))

    ## Change logical variables to TRUE/FALSE
    metadata$paired_end <- as.logical(metadata$paired_end)
    metadata$sra_misreported_paired_end <- as.logical(metadata$sra_misreported_paired_end)
    
    ## Load SRA-run info
    sra <- read.csv('/dcl01/leek/data/gtex_work/runs/sra/v2/hg38/SraRunInfo.csv',
        header = TRUE, stringsAsFactors = FALSE)
    colnames(sra) <- tolower(colnames(sra))

    ## Find average read length
    i <- match(metadata$run, sra$run)
    stopifnot(identical(metadata$run, sra$run[i]))
    metadata$avg_read_length <- sra$avglength[i]
} else if (opt$project == 'gtex') {
    stop('Have not implemented metadata for gtex project')
} else {
    stop("Invalid 'project' choice. Use gtex or sra")
}


## Find bigwig files
bigwigs <- system(paste0('cut -f 5 -d " " /dcl01/leek/data/gtex_work/runs/recount2/bwtool/bwtool_cmds_', opt$project, '.txt'), intern = TRUE)
names(bigwigs) <- gsub('.*coverage_bigwigs/|.bw', '', bigwigs)
j <- match(metadata$run, names(bigwigs))

## Matches number of bigwig files
stopifnot(sum(is.na(j)) == sum(is.na(metadata$auc)))
metadata$bigwig_path <- bigwigs[j]
metadata$bigwig_file <- gsub('.*coverage_bigwigs/', '', metadata$bigwig_path)

## Locate tsv files
tsv <- dir('/dcl01/leek/data/recount2/coverage', pattern = 'tsv', full.names = TRUE)
names(tsv) <- gsub('.sum.tsv', '', dir('/dcl01/leek/data/recount2/coverage', pattern = 'tsv'))
k <- match(metadata$run, names(tsv))

## Number of tsv files matches number of bigwig files
stopifnot(sum(is.na(k)) == sum(is.na(metadata$auc)))
metadata$tsv_path <- tsv[k]

## Save the metadata
save(metadata, file = paste0('metadata_', opt$project, '.Rdata'))

## Explore metadata
print('Number of dimensions')
dim(metadata)

print('Number of unique project IDs')
length(unique(metadata$project))

print('Number of unique run IDs')
length(unique(metadata$run))

print('First couple of rows')
head(metadata)

print('Number of NAs per column')
sapply(metadata, function(x) { sum(is.na(x)) })

print('Percent of NAs per column')
sapply(metadata, function(x) { sum(is.na(x)) }) / nrow(metadata) * 100

## Save a file per project ID
dir.create(paste0('project_metadata_', opt$project), showWarnings = FALSE)

metadata_clean <- metadata[, !colnames(metadata) %in% c('bigwig_path',
    'tsv_path')]
save(metadata_clean, file = paste0('metadata_clean_', opt$project, '.Rdata'))

meta <- split(metadata_clean, metadata_clean$project)
stopifnot(length(meta) == length(unique(metadata_clean$project)))

xx <- sapply(unique(metadata_clean$project), function(project) {
    project_metadata <- meta[[project]]
    write.table(project_metadata, file.path(paste0('project_metadata_',
        opt$project), paste0(project, '.tsv')), sep = '\t', row.names = FALSE,
        quote = FALSE)
})

## Save project ids in a file
write.table(unique(metadata_clean$project), file = paste0('project_ids_',
    opt$project, '.txt'), sep = '\t', row.names = FALSE, quote = FALSE,
    col.names = FALSE)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
