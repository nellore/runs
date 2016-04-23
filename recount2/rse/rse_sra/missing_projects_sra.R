## Usage
# module load R/3.3
# Rscript missing_projects_sra.R > missing_projects_sra_log.txt 2>&1

load('/dcl01/leek/data/gtex_work/runs/recount2/metadata/metadata_sra.Rdata')
missing <- c('ERP003617', 'ERP003731', 'ERP009290', 'ERP012951', 'SRP053791')
all(is.na(metadata$tsv_path[metadata$project %in% missing]))
missed <- metadata[metadata$project %in% missing, ]
nrow(missed)
missed

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
