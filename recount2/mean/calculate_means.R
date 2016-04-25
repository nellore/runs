## Load libraries
library('getopt')

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
dir.create(paste0('means_', opt$project), showWarnings = FALSE)

## Load metadata
load(opt$metadata)

## Subset to project of interest
metadata <- subset(metadata, project == opt$projectid)
if(nrow(metadata) == 0) stop(paste('Invalid project id', opt$projectid))

## Subset to only use the samples that have tsv files
metadata <- metadata[!is.na(metadata$tsv_path), ]
if(nrow(metadata) == 0) stop(paste('No samples have bwtool tsv files for project', opt$projectid))
rownames(metadata) <- NULL

## For testing
if(FALSE) {
    metadata <- metadata[1:3, ]
}

## Name resulting mean.bw file
outbw <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/mean', 
    paste0('means_', opt$project), paste0('mean_', opt$projectid, '.bw'))
outwig <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/mean', 
    paste0('means_', opt$project), paste0('mean_', opt$projectid, '.wig'))
    

scaleWig <- function(m) {
    paste(paste('scale', round(1e6*100*40 / m$auc, digits = 17),
        m$bigwig_path), collapse = ' ')
}
runCmd <- function(cmd, projectid, i = NULL) {
    if(is.null(i)) {
        shell_name <- paste0('.createWig_', opt$projectid, '.sh')
    } else {
        shell_name <- paste0('.createWig_', opt$projectid, '_part', i, '.sh')
    }    
    cat(cmd, file = shell_name)
    system(paste('sh', shell_name))
}

## Calculate mean bigwig
if(nrow(metadata) < 100) {
    ## Scale commands
    cmd <- scaleWig(metadata)
    ## Calculate mean wig file
    message(paste(Sys.time(), 'creating file', outwig))
    cmd <- paste('wiggletools write', outwig, 'mean', cmd)
    system.time( runCmd(cmd, opt$projectid) )
} else {
    ## Define subsets to work on
    sets <- Hmisc::cut2(seq_len(nrow(metadata)), m = 50)
    meta <- split(metadata, sets)
    names(meta) <- seq_len(length(meta))
    
    ## Calculate sums per subsets
    system.time( tmpfiles <- mapply(function(m, i) {
        cmd <- scaleWig(m)
        
        ## Use TMPDIR if available
        if(Sys.getenv('TMPDIR') == '') {
            tmpdir <- file.path('/dcl01/leek/data/gtex_work/runs/recount2/mean', paste0('means_', opt$project), paste0('.', opt$projectid))
            dir.create(tmpdir, showWarnings = FALSE)
        } else {
            tmpdir <- tempdir()
        }        
        tmpwig <- file.path(tmpdir, paste0('sum_', opt$projectid, '_part',
            i, '.wig'))
        message(paste(Sys.time(), 'creating file', tmpwig))
        cmd <- paste('wiggletools write', tmpwig, 'sum', cmd)
        runCmd(cmd, opt$projectid, i)
        return(tmpwig)
    }, meta, names(meta)) )
    
    ## Calculate final mean
    cmd <- paste('wiggletools write', outwig, 'scale', 1/nrow(metadata), 'sum', 
        paste(tmpfiles, collapse = ' '))
    system.time( runCmd(cmd, opt$projectid) )
}

## Transform to bigwig file
message(paste(Sys.time(), 'creating file', outbw))
cmd2 <- paste('wigToBigWig', outwig,
    '/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', outbw)
system.time( system(cmd2) )

## Testing: compare with R
if(FALSE) {
    library('rtracklayer')
    library('GenomicRanges')
    library('derfinder')
    bws <- BigWigFileList(metadata$bigwig_path)
    chrinfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes',
        header = FALSE, col.names = c('chr', 'length'), sep = '\t',
        stringsAsFactors = FALSE)
    #chrs <- chrinfo$chr # not all of them present in all bigwigs
    chrs <- paste0('chr', c(1:22, 'X', 'Y'))
    chrlen <- chrinfo$chrs[chrinfo$chr %in% chrs]
    fullCov <- fullCoverage(bws, chrs, chrlens = chrlen,
        totalMapped = metadata$auc, targetSize = 1e6*100*40)
    
    mean <- lapply(fullCov, function(x) { Reduce('+', x) / ncol(x)})
    
    meanBW <- import(outbw, as = 'RleList')
    
    differences <- lapply(chrs, function(chr) {
        mean[[chr]] - meanBW[[chr]]
    })
    names(differences) <- chrs
    
    ## Are numerical differences small?
    num_small <- lapply(differences, function(x) { abs(x) < 0.0153 })
    num_small
    all(unlist(sapply(num_small, runValue)))
    
    ## Summary of differences
    diff_summary <- lapply(differences, summary)
    diff_summary
}


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
