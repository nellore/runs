##
library(Matrix)
library(parallel)
library(GenomicRanges)

# ######### phenotype ###################

# read in phenotype data, first drop column 89
thecall = pipe("cut -f1-88,90- /dcl01/leek/data/sraintrons/all_SRA_metadata.tsv")
pd = read.delim(thecall, as.is=TRUE)
rownames(pd) = pd$run_accession

## add abstract column
absCall = pipe("cut -f89 /dcl01/leek/data/sraintrons/all_SRA_metadata.tsv")
ab = scan(absCall, sep="\n", what="character")
pd$sra_study_abstract = ab[-1]

##################################
##### designate count files ######
chrFiles = paste0("/dcl01/lieber/ajaffe/AllSRA/byChr/chr",
	c(1:22,"X","Y", "M"), "_junctions.rda")
names(chrFiles) = paste0("chr", c(1:22,"X","Y", "M"))

# read in counts and filter
N = 1000 # number of samples with non-zero expression
countsFiltered = lapply(chrFiles, function(fn) {
	cat(".")
	# load data for that chromosome
	vars = load(fn)
	
	# rename 
	newvars = c("counts", "map", "adjCounts")
	for(i in seq(along=vars)) assign(newvars[i], get(vars[i]))
	
	# mean adjusted coverage more than 0
	exprsIndex = which(rowSums(counts > 0) >= N)
	
	# filter
	mapFilter = map[exprsIndex]
	countsFilter = counts[exprsIndex,]

	outList = list(counts = countsFilter, map = mapFilter)
	return(outList)
})

####### match up filtered counts
countList = lapply(countsFiltered, function(x) x$counts[,rownames(pd)])
mapList = GRangesList(lapply(countsFiltered, function(x) x$map))

sum(sapply(countList, nrow)) # 305k

save(pd, countList, mapList, file="rdas/filtered_sra_data_N1000_rawCounts.rda")
