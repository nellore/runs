library('GenomicRanges')

message(paste(Sys.time(), 'reading regions -- top 500 files'))
system.time(load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/regions-cut0.5.Rdata'))
print(object.size(regions), units = 'Mb')
message(paste(Sys.time(), 'reading coverage matrix -- top 500 files'))
system.time(load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/coverageMatrix-cut0.5.Rdata'))
print(object.size(coverageMatrix), units = 'Gb')

regs <- regions
covMat <- coverageMatrix


message(paste(Sys.time(), 'reading regions -- all files'))
system.time(load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/bwtool/regions-cut0.5.Rdata'))
print(object.size(regions), units = 'Mb')
message(paste(Sys.time(), 'reading coverage matrix -- all files'))
system.time(load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/bwtool/coverageMatrix-cut0.5.Rdata'))
print(object.size(coverageMatrix), units = 'Gb')

mcols(regs) <- NULL
identical(regs, regions)

identical(colnames(covMat), colnames(coverageMatrix)[1:500])
identical(covMat, coverageMatrix[, 1:500])
summary(covMat - coverageMatrix[, 1:500])

diffs <- lapply(seq_len(10), function(i) {
    summary(covMat[, i] - coverageMatrix[, i])
})

which_diff <- lapply(seq_len(10), function(i) {
    which(covMat[, i] != coverageMatrix[, i])
})
summary(sapply(which_diff, length))
summary(sapply(which_diff, length)) / nrow(coverageMatrix) * 100

nas <- sapply(seq_len(10), function(x) { which(is.na(coverageMatrix[, x]))})
regions[nas[[4]]]
colnames(coverageMatrix)[4]

## Explore differences in first sample
covMat[head(which_diff[[1]]), 1]
coverageMatrix[head(which_diff[[1]]), 1]

## It's not a constant difference
covMat[head(which_diff[[1]]), 1] / coverageMatrix[head(which_diff[[1]]), 1]

## Is it region length?
covMat[head(which_diff[[1]]), 1]
coverageMatrix[head(which_diff[[1]]), 1] * width(regions[head(which_diff[[1]])])

covMat[head(which_diff[[1]]), 1] / (coverageMatrix[head(which_diff[[1]]), 1] * width(regions[head(which_diff[[1]])]))

coverageMatrix[head(which_diff[[1]]), 1] * (width(regions[head(which_diff[[1]])]))


table(row.names(covMat))

summary(covMat[, 1] - coverageMatrix[, 1])
summary(sort(covMat[, 1]) - sort(coverageMatrix[, 1]))

regions[which_diff[[1]]]


## Get the raw counts for regions 1, 2, 6, 7, 9, 12 from sample 1 which
## is /dcl01/leek/data/gtex/batch_19/coverage_bigwigs/SRR660824_SRS389722_SRX222703_male_lung.bw
library('rtracklayer')
diff_1 <- head(which_diff[[1]])
regs_1 <- regions[diff_1]

export(regs_1, con  = '/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/misc/regions_debug.BED', format = 'BED')

load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/pheno/pheno_missing_less_10.Rdata')
sampleFile <- pheno$BigWigPath[1]
info <- import(sampleFile, selection = regs_1, as = 'RleList')[unique(seqnames(regs_1))]

sums_1 <- sapply(as.character(unique(seqnames(regs_1))), function(chr) {
    sum(Views(info[[chr]], ranges(regs_1[seqnames(regs_1) == chr])))
})
unlist(sums_1)



## Find count files
counts_files <- file.path(dir('/dcl01/leek/data/gtex', pattern = 'batch', full.names = TRUE), 'cross_sample_results', 'counts.tsv.gz')
names(counts_files) <- dir('/dcl01/leek/data/gtex', pattern = 'batch')

## Read in counts info
counts <- lapply(counts_files, function(i) {
    read.table(i, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
})
counts <- do.call(rbind, counts)
counts$totalMapped <- as.integer(sapply(strsplit(counts$total.mapped.reads, ','), '[[', 1))

## Match files to counts
map <- match(gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', sampleFile), counts$X)
counts <- counts[map, ]

## Compare manual calculation done with R with railMatrix() R results
covMat[diff_1, 1] - unlist(sums_1) / (counts$totalMapped / 40e6) / (pheno$avgLength[1] / 2)


unlist(sums_1)
counts$totalMapped
pheno$avgLength[1]
unlist(sums_1) / (counts$totalMapped / 40e6) / (pheno$avgLength[1] / 2)


sampleNames <- gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', sampleFile)

## Read in bwtool results for sample 1
## (bwtool using sum, total # mapped not total # reads)
cov_sample1 <- mapply(function(tsvFile, sampleName) {
    message(paste(Sys.time(), 'reading file', tsvFile))
    res <- read.table(tsvFile, header = FALSE, colClasses = list(NULL, NULL, NULL, 'numeric'))
    colnames(res) <- sampleName
    return(as.matrix(res))
}, '/dcl01/leek/data/gtex_work/cov_sums/SRR660824.sum.tsv', sampleNames, SIMPLIFY = FALSE)
cov_sample1 <- do.call(cbind, cov_sample1)

## Explore differences between sample 1  and R results
identical(cov_sample1, covMat[, 1])
identical(colnames(cov_sample1), colnames(covMat)[1])
differences <- cov_sample1[, 1] - covMat[, 1]
summary(differences)
diff_i <- which(abs(differences) > 1e-4)
length(diff_i)
length(which(abs(differences) > 1e-8))
length(which(abs(differences) > 1e-9))
differences[head(diff_i)]
regs[head(diff_i)]
cov_sample1[head(diff_i), 1]
covMat[head(diff_i), 1]

## Explore the larger differnces
length(which(abs(differences) > 0.1))
differences[which(abs(differences) > 0.1)]
cov_sample1[which(abs(differences) > 0.1), 1]
covMat[which(abs(differences) > 0.1), 1]
regs[which(abs(differences) > 0.1)]

## Calculate manually coverage for these 9 regions
regs_2 <- regs[which(abs(differences) > 0.1)]
info_2 <- import(sampleFile, selection = regs_2, as = 'RleList')[unique(seqnames(regs_2))]

sums_2 <- sapply(as.character(unique(seqnames(regs_2))), function(chr) {
    sum(Views(info_2[[chr]], ranges(regs_2[seqnames(regs_2) == chr])))
})
unlist(sums_2)

## Numbers used
counts$totalMapped
pheno$avgLength[1]
unlist(sums_2) / (counts$totalMapped / 40e6) / (pheno$avgLength[1] / 2)

(unlist(sums_2) / (counts$totalMapped / 40e6) / (pheno$avgLength[1] / 2))[5]

res_2 <- covMat[which(abs(differences) > 0.1), 1] - unlist(sums_2) / (counts$totalMapped / 40e6) / (pheno$avgLength[1] / 2)
res_2
res_2[5]

