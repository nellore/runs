# mkdir -p logs
# Rscript create_gtexChr.R > logs/create_gtexChr_log.txt 2>&1
library('GenomicRanges')
library('devtools')

## Load chr matching info
chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/genomicState/hg38.ucsc.sizes.ensembl.gencode', header = TRUE, colClasses = c('character', 'numeric', 'character', 'character'))


## Process memory used
mergeLog <- readLines('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/bwtool/logs/gtex-merge-bwtool.e8877097')

mem <- mergeLog[grep('coverage matrix for.*uses', mergeLog)]

memInfo <- as.numeric(gsub('.*uses | Mb in RAM| Gb in RAM', '', mem))
names(memInfo) <- gsub('.*for | uses.*', '', mem)
mem[grep('Gb in RAM', mem)]

m <- match(chrInfo$ucsc, names(memInfo))
m_clean <- m[!is.na(m)]
chrInfo$RAM <- NA
chrInfo$RAM[!is.na(m)] <- memInfo[m_clean]
chrInfo$RAM_unit <- NA
chrInfo$RAM_unit[which(!is.na(m))[grep('Gb in RAM', mem)]] <- 'Gb'
chrInfo$RAM_unit[which(!is.na(m))[grep('Mb in RAM', mem)]] <- 'Mb'

chrInfo$matrixFile[!is.na(m)] <- file.path('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/bwtool', paste0('coverageMatrix-cut0.5-', chrInfo$ucsc[!is.na(m)], '.Rdata'))

load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/regions-cut0.5.Rdata')

n_regs <- sapply(split(regions, seqnames(regions)), length)

chrInfo$nRegions <- NA
chrInfo$nRegions[!is.na(m)] <- n_regs[match(chrInfo$ucsc[!is.na(m)], names(n_regs))]

write.table(chrInfo, file = 'gtexChr.txt', sep = '\t', row.names = FALSE, quote = FALSE)

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
