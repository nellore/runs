## Prepare metadata
# module load R/3.3
# mkdir -p logs
# Rscript create_meta_web.R -p "sra" > logs/create_meta_web_sra_log.txt 2>&1
# Rscript create_meta_web.R -p "gtex" > logs/create_meta_web_gtex_log.txt 2>&1

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

## For testing
if(FALSE) {
    opt <- list(project = 'sra')
}

## Load metadata
file <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/metadata/',
    paste0('metadata_', opt$project, '.Rdata'))
stopifnot(file.exists(file))
load(file)

## Load abstracts
file <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/metadata_abstract/', 
    paste0('abstracts_', opt$project, '.Rdata'))
stopifnot(file.exists(file))
load(file)

projects <- unique(metadata$project)
meta_web <- data.frame(
    accession = paste0(
        '<a href="http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=',
        projects, '">', projects, '</a>'),
    number_samples = sapply(projects, function(project) {
        sum(metadata$project == project)}),
    species = 'human',
    abstract = abstracts$study_abstract[match(projects,
        abstracts$study_accession)],
    rse_gene = NA,
    rse_exon = NA,
    counts_gene = NA,
    counts_exon = NA,
    phenotype = NA,
    genes = '<a href="https://lcolladotor.shinyapps.io/recount/ucsc-knowngene-hg38-genes-bp-length.Rdata">link</a>',
    exons = '<a href="https://lcolladotor.shinyapps.io/recount/ucsc-knowngene-hg38-exons.Rdata">link</a>',
    files_info = NA,
    stringsAsFactors = FALSE
)
rownames(meta_web) <- NULL
for(project in projects) {
    if(dir.exists(file.path('/dcl01/leek/data/gtex_work/runs/recount2/rse/', paste0('rse_', opt$project), project))) {
        ## Have to change this to actual URLs once the data is uploaded        
        meta_web$rse_gene[projects == project] <- paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/rse_gene.Rdata">link</a>')
        meta_web$rse_exon[projects == project] <- paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/rse_exon.Rdata">link</a>')
        meta_web$counts_gene[projects == project] <- paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/counts_gene.tsv.gz">link</a>')
        meta_web$counts_exon[projects == project] <- paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/counts_exon.tsv.gz">link</a>')
    }
    if(file.exists(file.path('/dcl01/leek/data/gtex_work/runs/recount2/metadata/', paste0('project_metadata_', opt$project), paste0(project, '.tsv')))) {
       ## Have to change this to actual URLs once the data is uploaded    
        meta_web$phenotype[projects == project] <- paste0(
            '<a href="http://duffel.rail.bio/recount/', project, '/', project,
            '.tsv">link</a>')
    }
    if(dir.exists(file.path('/dcl01/leek/data/gtex_work/runs/recount2/fileinfo/', paste0('fileinfo_', opt$project), project))) {
        ## Might host these files on the website itself, if so we'll have to
        ## change them to relative URLs
        meta_web$files_info[projects == project] <- paste0(
            '<a href="http://duffel.rail.bio/recount/', project,
            '/files_info.tsv">link</a>')
    }
}

## Explore number of samples per project
message(paste(Sys.time(), "number of samples per project"))
summary(meta_web$number_samples)

## How much of it is missing?
message(paste(Sys.time(), "number missing"))
sapply(meta_web, function(x) sum(is.na(x)))

## Save info
save(meta_web, file = paste0('meta_web_', opt$project, '.Rdata'))
write.table(meta_web, file = paste0('meta_web_', opt$project, '.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
