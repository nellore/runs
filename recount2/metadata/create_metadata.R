## Prepare metadata
# module load R/3.3
# mkdir -p logs
# Rscript create_metadata.R > logs/create_metadata.txt 2>&1

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

## Find bigwig files
bigwigs <- system('cut -f 5 -d " " /dcl01/leek/data/gtex_work/runs/recount2/bwtool/bwtool_cmds_sra.txt', intern = TRUE)
names(bigwigs) <- gsub('.*coverage_bigwigs/|.bw', '', bigwigs)
j <- match(metadata$run, names(bigwigs))

## Matches number of bigwig files
stopifnot(sum(is.na(j)) == sum(is.na(metadata$auc)))
metadata$bigwig_path <- bigwigs[j]
metadata$bigwig_file <- gsub('.*coverage_bigwigs/', '', metadata$bigwig_path)

## Save the metadata
save(metadata, file = 'metadata.Rdata')

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
dir.create('project_metadata', showWarnings = FALSE)

metadata_nobw <- metadata[, colnames(metadata) != 'bigwig_path']
save(metadata_nobw, file = 'metadata_nobw.Rdata')

meta <- split(metadata_nobw, metadata_nobw$project)
stopifnot(length(meta) == length(unique(metadata_nobw$project)))

xx <- sapply(unique(metadata_nobw$project), function(project) {
    project_metadata <- meta[[project]]
    write.table(project_metadata, file.path('project_metadata', paste0(project, '.tsv')), sep = '\t', row.names = FALSE, quote = FALSE)
})

## Save project ids in a file
write.table(unique(metadata_nobw$project), file = 'project_ids.txt', sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
