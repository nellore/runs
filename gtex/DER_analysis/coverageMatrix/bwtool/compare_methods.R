library('GenomicRanges')
system.time(load('regions-cut0.5.Rdata'))
system.time(load('coverageMatrix-cut0.5.Rdata'))
print(object.size(coverageMatrix), units = 'Gb')

regs <- regions
covMat <- coverageMatrix

system.time(load('bwtool/regions-cut0.5.Rdata'))
system.time(load('bwtool/coverageMatrix-cut0.5.Rdata'))

mcols(regs) <- NULL
identical(regs, regions)

identical(colnames(covMat), colnames(coverageMatrix)[1:500])
identical(covMat, coverageMatrix[, 1:500])
summary(covMat - coverageMatrix[, 1:500])

diffs <- lapply(seq_len(500), function(i) {
    summary(covMat[, i] - coverageMatrix[, i])
})