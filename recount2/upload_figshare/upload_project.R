## Load libraries
library('getopt')
library('rfigshare')
library('tools')
library('RCurl')

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
    opt <- list(project = 'sra', projectid = 'DRP000499',
        'metadata' = '/dcl01/leek/data/gtex_work/runs/recount2/metadata/metadata_sra.Rdata'
    )
    ## Largest one, to find memory needed
    opt <- list(project = 'sra', projectid = 'SRP025982',
        'metadata' = '/dcl01/leek/data/gtex_work/runs/recount2/metadata/metadata_sra.Rdata'
    )
}

## Create output dir
dir.create(paste0('upload_', opt$project), showWarnings = FALSE)

## Load metadata
load(opt$metadata)

## Subset to project of interest
metadata <- subset(metadata, project == opt$projectid)
if(nrow(metadata) == 0) stop(paste('Invalid project id', opt$projectid))

## Subset to only use the samples that have tsv files
metadata <- metadata[!is.na(metadata$tsv_path), ]
if(nrow(metadata) == 0) stop(paste('No samples have bwtool tsv files for project', opt$projectid))
rownames(metadata) <- NULL

## Locate counts and rse files
rse_path <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/rse',
    paste0('rse_', opt$project), opt$projectid)
rse_files <- dir(rse_path, full.names = TRUE)
names(rse_files) <- dir(rse_path)

rse_up <- c('counts_exon.tsv.gz', 'counts_gene.tsv.gz', 'rse_exon.Rdata',
    'rse_gene.Rdata')
if(!all(rse_up %in% names(rse_files))) stop(paste('Missing counts/rse files for project',
    opt$projectid))

## Locate metadata file
meta_path <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/metadata',
    paste0('project_metadata_', opt$project), paste0(opt$projectid, '.tsv'))
names(meta_path) <- paste0(opt$projectid, '.tsv')
if(!file.exists(meta_path)) stop(paste('Missing metadata file for project', opt$projecid))

## Find mean bigwig file
mean_bw <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/mean', 
    paste0('means_', opt$project), paste0('mean_', opt$projectid, '.bw'))
if(!file.exists(mean_bw)) stop(paste('Missing mean bigwig file for project', opt$projecid))
names(mean_bw) <- paste0('mean_', opt$projectid, '.bw')

## Create output dir for the project
outdir <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/upload_figshare', 
    paste0('upload_', opt$project), opt$projectid)
dir.create(outdir, showWarnings = FALSE)
stopifnot(dir.exists(outdir))

## Find bigwig files
upload_bigwig <- metadata$bigwig_path
names(upload_bigwig) <- metadata$bigwig_file


## List of files to upload
upload_files <- c(rse_files[rse_up], meta_path, mean_bw, upload_bigwig)

## Compute md5sums
print('Time spent calculating md5sums')
system.time( files_md5sum <- md5sum(upload_files) )

## Save md5sums
write.table(data.frame(file = names(upload_files), md5sum = files_md5sum),
    file = file.path(outdir, 'files_md5sum.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)

## Final list of files to uplaod
upload_files <- c(upload_files, 'files_md5sum.tsv' = file.path(outdir,
    'files_md5sum.tsv'))

## Calculate total file size
print('Time spent calculating file sizes')
system.time( file_size <- as.numeric(system(paste('du -l',
    paste(upload_files, collapse = ' '), '| cut -f1'), intern = TRUE)) )
names(file_size) <- names(upload_files)
project_size <- sum(file_size / 1024)
if(project_size < 1024) {
    print(paste('The total file size for project', opt$projectid, 'is', round(project_size, digits = 3), 'Mb.'))
} else {
    print(paste('The total file size for project', opt$projectid, 'is', round(project_size / 1024, digits = 3), 'Gb.'))
}

## Create article
art_id <- fs_create(
    title = paste('recount2: files for project', opt$projectid),
    description = paste0(
        'This figshare fileset is part of the recount2 project. It contains bigWig files for samples aligned with Rail-RNA for project ', opt$projectid, ' as well as RangedSummarizedExperiment R objects for the gene-level and exon-level raw counts (files rse_gene.Rdata and rse_exon.Rdata). The counts and metadata used for the RangedSummarizedExperiment objects is also included in this fileset (files counts_exons.tsv.gz, counts_exons.tsv.gz and ', opt$projectid, '.tsv). The file mean_', opt$project_id, '.bw can be used with http://bioconductor.org/packages/derfinder to perform a differential expression analysis at the expressed regions-level. Finally, files_md5sum.tsv and files_url.tsv contain the md5sum hashes and download URLs for the files, respectively. You can find other projects via the recount2 website at ADD_URL. Using http://bioconductor.org/packages/recount you can download the data. The code for processing these files can be found at https://github.com/nellore/runs. If you use this fileset please cite the recount paper ADD_DOI and the DOI for this fileset hosted at figshare. Thank you!'),
    type = 'fileset'
)
write.table(art_id, file = file.path(outdir, 'article_id.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
    
## Add category, recount2 tag, and link to recount2 website
fs_add_recount2 <- function (article_id) {
    o <- fs_details(article_id, mine = TRUE)
    tag <- "recount2 project"
    m <- match(tag, sapply(o$tags, `[[`, "name"))
    out <- if (is.na(m)) 
        fs_add_tags(article_id, tag)
    invisible(out)
}
fs_add_categories(art_id, 'Computational Biology')
fs_add_links(art_id, c('http://RECOUNT2.org',
    'http://bioconductor.org/packages/recount',
    'https://github.com/nellore/runs',
    'http://rail.bio/',
    'http://bioconductor.org/packages/derfinder')
)
fs_add_recount2(art_id)

## If the project is too large, make it public first?
#if(project_size > )

## Load token info
load('.auth_string.Rdata') ## String from creating a "Personal Token" at https://figshare.com/account/applications

## Upload files
fs_my_upload <- function(article, files, token) {
    cmd <- paste('python /dcl01/leek/data/gtex_work/runs/recount2/figshare.py --article-id', article, '--paths', paste(files, collapse = ' '), '--token', token)
    system(cmd)
}
print('Time spent uploading files')
system.time( upload_info <- fs_my_upload(art_id, upload_files,
    token = auth_string) )
#save(upload_info, file = file.path(outdir, 'upload_info.Rdata'))

## Make public
fs_make_public(art_id)

## Get short url #http://stackoverflow.com/questions/6500721/find-where-a-t-co-link-goes-to
unshorten_url <- function(uri){
    if(RCurl::url.exists(uri)){
        # listCurlOptions()
        opts <- list(
            followlocation = TRUE,  # resolve redirects
            ssl.verifyhost = FALSE, # suppress certain SSL errors
            ssl.verifypeer = FALSE, 
            nobody = TRUE, # perform HEAD request
            verbose = FALSE
        )
        curlhandle <- getCurlHandle(.opts = opts)
        getURL(uri, curl = curlhandle)
        info <- getCurlInfo(curlhandle)
        rm(curlhandle)  # release the curlhandle!
        info$effective.url
    } else {
        # just return the url as-is
        uri
    }
}

## Get list of files and URLs
details <- fs_details(art_id, mine = FALSE)
down_info <- fs_download(art_id)
url_info <- data.frame(url = down_info[, 1], 
    name = sapply(details$files, '[[', 'name'),
    stringsAsFactors = FALSE
)
url_info$true_url <- sapply(url_info$url, unshorten_url)
write.table(url_info,
    file = file.path(outdir, 'files_url.tsv'),
    sep = '\t', row.names = FALSE, quote = FALSE, col.names = TRUE)

## Upload files info
fs_upload(art_id, file.path(outdir, 'files_url.tsv'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
