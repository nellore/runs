######
library(Matrix)
library(parallel)
library(GenomicRanges)

## combine into 1 matrix?
load("rdas/filtered_sra_data_N1000_rawCounts.rda")
jMap = unlist(mapList,use.names=FALSE)

## library size normalize
pd$totalReads = pd$sra_spots
pairIndex=which(pd$sra_layout == "PAIRED")
pd$totalReads[pairIndex] = 2*pd$sra_spots[pairIndex]

pd$sumFilterJxnCount = rowSums(sapply(countList, colSums))
f = lm(totalReads ~ sumFilterJxnCount, data=pd,
	subset=totalReads > 10000)
pd$totalReadsPred = pd$totalReads
pd$totalReadsPred[pd$totalReads < 1000] = predict(f, 
	pd[pd$totalReads < 1000,])

## make matrix
library(parallel)
matList = lapply(countList, as.matrix)	
pd$mappedPer10M = pd$totalReadsPred/10e6
matListAdj = lapply(matList, function(x) {
	cat(".")
	cc = matrix(pd$mappedPer10M, 
		nc = ncol(x), nr = nrow(x),byrow=TRUE)
	x/cc
})
jMap$meanExprs = unlist(lapply(matList, rowMeans))

matListAdj = lapply(matListAdj, function(x) {
	cat(".")
	log2(x+1)
})

adjCounts = do.call("rbind", matListAdj)
identical(rownames(adjCounts), names(jMap)) # TRUE

save(adjCounts, jMap, pd,
	file="rdas/filtered_sra_data_N1000_normCounts.rda")
