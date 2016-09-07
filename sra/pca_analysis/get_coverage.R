library(readr)
library(Matrix)
library(pbapply)
library(GenomicRanges)
library(pbapply)
load("/dcl01/leek/data/sraintrons/filtered_sra_data_N1000_rawCounts.rda")

sumsList <- pblapply(countList, colSums, na.rm=TRUE)
sumsList <- do.call(rbind, sumsList)
sumsList <- colSums(sumsList)
totalJunctionCovRaw <- sumsList

save(totalJunctionCovRaw, file="extdata/coverage_raw_info.rda")