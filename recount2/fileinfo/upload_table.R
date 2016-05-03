## Usage
# mkdir -p logs
# module load R/3.3
# Rscript upload_table.R > logs/upload_table_log.txt 2>&1

## Identify files to upload
upload <- dir('/dcl01/leek/data/gtex_work/runs/recount2/fileinfo/fileinfo_sra', full.names = TRUE)
upload <- file.path(upload, 'upload_files.Rdata')
names(upload) <-  dir('/dcl01/leek/data/gtex_work/runs/recount2/fileinfo/fileinfo_sra')

## Find all the info to upload
upload_table <- mapply(function(rda, project) {
    load(rda)
    data.frame(path = upload_files, file_name = names(upload_files),
        project = project, stringsAsFactors = FALSE)
}, upload, names(upload), SIMPLIFY = FALSE)
upload_table <- do.call(rbind, upload_table)
rownames(upload_table) <- NULL

## Explore table a little bit
head(upload_table)
dim(upload_table)

## Save info
save(upload_table, file = 'upload_table.Rdata')

write.table(upload_table, file = 'upload_table.tsv', sep = '\t',
    row.names = FALSE, quote = FALSE, col.names = FALSE)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
