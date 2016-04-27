## Prepare metadata
# module load R/3.3
# mkdir -p logs
# Rscript get_abstract.R -p "sra" > logs/get_abstract_sra_log.txt 2>&1
# Rscript get_abstract.R -p "gtex" > logs/get_abstract_gtex_log.txt 2>&1

library('getopt')
library('SRAdb')

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

## For testing
if(FALSE) {
    opt <- list(project = 'sra')
}

## Load metadata
file <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/metadata/',
    paste0('metadata_', opt$project, '.Rdata'))
stopifnot(file.exists(file))
load(file)

sqlfile <- 'SRAmetadb.sqlite'
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()

## Create connection
sra_con <- dbConnect(SQLite(), sqlfile)

## Show file info
message(paste(Sys.time(), "information about", sqlfile))
file.info(sqlfile)
dat <- dbGetQuery(sra_con, "SELECT * FROM metaInfo")
dat

## Get unique projects
projects <- unique(metadata$project)

## Find the abstract for each project
abstracts <- dbGetQuery(sra_con,
    paste0("SELECT study_abstract, study_accession FROM study WHERE study_accession IN ('",
    paste(projects, collapse="', '"), "')")
)

## Number missing
message(paste(Sys.time(), "number missing"))
sapply(abstracts, function(x) sum(is.na(x)))

## Save info
save(abstracts, file = paste0('abstracts_', opt$project, '.Rdata'))

## Close connection
dbDisconnect(sra_con)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
