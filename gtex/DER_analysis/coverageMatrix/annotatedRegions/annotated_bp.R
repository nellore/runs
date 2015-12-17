# mkdir -p logs
# Rscript annotated_bp.R > logs/annotated_bp_log.txt 2>&1
library('GenomicRanges')
library('devtools')

load('annotated_ensembl.Rdata')
load('annotated_ensembl_8.Rdata')
load('annotated_ensembl_20.Rdata')
load('annotated_ucsc.Rdata')
load('annotated_ucsc_8.Rdata')
load('annotated_ucsc_20.Rdata')
load('regs_ensembl.Rdata')
regs <- regs_ensembl

annotated <- list(
    'ensembl' = list('1bp' = annotated_ensembl, '8bp' = annotated_ensembl_8, '20bp' = annotated_ensembl_20),
    'ucsc' = list('1bp' = annotated_ucsc, '8bp' = annotated_ucsc_8, '20bp' = annotated_ucsc_20)
)

annotation_bp <- lapply(names(annotated), function(db) {
    resov <- lapply(c('1bp', '8bp', '20bp'), function(minov) {
        reslen <- lapply(c(0, 7, 19), function(minlen) {
            s <- width(regs) > minlen
            res <- sapply(c('exon', 'intergenic', 'intron'), function(type, subsetIndex = s) {
                sum(width(regs)[annotated[[db]][[minov]]$countTable[[type]] > 0 & s])
            })
            data.frame(t(res), 'minlen' = minlen, 'minov' = minov, db = db, stringsAsFactors = FALSE)
        })
        do.call(rbind, reslen)
    })
    do.call(rbind, resov)
})
annotation_bp <- do.call(rbind, annotation_bp)
annotation_bp$total_bp <- annotation_bp$exon + annotation_bp$intergenic + annotation_bp$intron
annotation_bp$exon_per <- annotation_bp$exon / annotation_bp$total_bp * 100
annotation_bp$intergenic_per <- annotation_bp$intergenic / annotation_bp$total_bp * 100
annotation_bp$intron_per <- annotation_bp$intron / annotation_bp$total_bp * 100
options(width = 120)
print(annotation_bp, digits = 4)


## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
