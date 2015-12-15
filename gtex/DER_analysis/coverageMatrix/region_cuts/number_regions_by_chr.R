library('GenomicRanges')
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/region_cuts/region_cuts.Rdata')
s <- split(region_cuts[['0.5']], seqnames(region_cuts[['0.5']]))
elementLengths(s)
