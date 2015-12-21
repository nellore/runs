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

cases <- list(
    'exon' = c(1, 0, 0),
    'intergenic' = c(0, 1, 0),
    'intron' = c(0, 0, 1),
    'exon_intron' = c(1, 0, 1),
    'exon_intergenic' = c(1, 1, 0),
    'intron_intergenic' = c(0, 1, 1),
    'all' = c(1, 1, 1),
    'none' = c(0, 0, 0)
)

annotation_bp <- lapply(names(annotated), function(db) {
    resov <- lapply(c('1bp', '8bp', '20bp'), function(minov) {
        reslen <- lapply(c(0, 7, 19), function(minlen) {
            s <- width(regs) > minlen
            res <- sapply(names(cases), function(type, subsetIndex = s) {
                count <- annotated[[db]][[minov]]$countTable[s, ] > 0
                case <- matrix(rep(cases[[type]] > 0, each = nrow(count)), ncol = 3)
                case_true <- apply(count == case, 1, all)
                sum(width(regs)[which(s)[case_true]])
            })
            data.frame(t(res), 'minlen' = minlen, 'minov' = minov, db = db, stringsAsFactors = FALSE)
        })
        do.call(rbind, reslen)
    })
    do.call(rbind, resov)
})
annotation_bp <- do.call(rbind, annotation_bp)
annotation_bp$total_bp <- rowSums(annotation_bp[, 1:8])
options(width = 160)
annotation_bp

annotation_bp$exon_per <- annotation_bp$exon / annotation_bp$total_bp * 100
annotation_bp$intergenic_per <- annotation_bp$intergenic / annotation_bp$total_bp * 100
annotation_bp$intron_per <- annotation_bp$intron / annotation_bp$total_bp * 100
annotation_bp$exon_intron_per <- annotation_bp$exon_intron / annotation_bp$total_bp * 100
annotation_bp$exon_intergenic_per <- annotation_bp$exon_intergenic / annotation_bp$total_bp * 100
annotation_bp$intron_intergenic_per <- annotation_bp$intron_intergenic / annotation_bp$total_bp * 100
annotation_bp$all_per <- annotation_bp$all / annotation_bp$total_bp * 100
annotation_bp$none_per <- annotation_bp$none / annotation_bp$total_bp * 100

print(annotation_bp[, 9:ncol(annotation_bp)], digits = 4)

## Save results
save(annotation_bp, file = 'annotation_bp.Rdata')

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
