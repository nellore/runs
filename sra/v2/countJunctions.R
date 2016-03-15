# qsub -V -l mf=40,h_vmem=60G,h_fsize=50G,h_stack=256M -t 1-25 -cwd -b y R CMD BATCH --no-save countJunctions.R
# qsub -V -l mf=200G,h_vmem=250G,h_fsize=50G,h_stack=256M -cwd -b y R CMD BATCH --no-save countJunctions.R

library(GenomicRanges)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg19)

### gunzip /dcl01/leek/data/sraintrons/all_SRA_introns.tsv.gz
fn = '/dcl01/lieber/ajaffe/AllSRA/all_SRA_introns.tsv'

## make annotation
# chr = read.delim(pipe(paste("cut -f1", fn)),header=FALSE)$V1
# start = read.delim(pipe(paste("cut -f2", fn)),header=FALSE)$V1
# end = read.delim(pipe(paste("cut -f3", fn)),header=FALSE)$V1
# strand = read.delim(pipe(paste("cut -f4", fn)),header=FALSE)$V1
# leftMotif = read.delim(pipe(paste("cut -f5", fn)),header=FALSE)$V1
# rightMotif = read.delim(pipe(paste("cut -f6", fn)),header=FALSE)$V1
# motif = paste0(leftMotif, ":", rightMotif)

# jMap = GRanges(chr, IRanges(start, end), strand=strand,
	# motif = motif)
# names(jMap) = paste0(chr, ":", start, "-", end, "(",strand, ")")
# seqlevels(jMap)= seqlevels(Hsapiens)[1:25]
# seqlengths(jMap) = seqlengths(Hsapiens)[1:25]
# genome(jMap) = "hg19"

# save(jMap, file = "junctionMap_allSRA.rda")

### read in annotation junction map
load("junctionMap_allSRA_annotated.rda")

# split up by chr
chrIndexes = split(seq(along=jMap), seqnames(jMap))
rowIndex= as.data.frame(t(sapply(chrIndexes, range)))
colnames(rowIndex) = c("rowStart", "rowEnd")
rowIndex$skip = rowIndex$rowStart-1
rowIndex$numRow = rowIndex$rowEnd - rowIndex$rowStart+1

## read in pheno data
pheno = read.delim("/dcl01/leek/data/sraintrons/index_to_SRA_accession.tsv",
	as.is=TRUE, header=FALSE)
colnames(pheno) = c("colIndex", "sraProject", 
	"sraSample", "sraExperiment", "sraRun")

## merge in total mapped
cc = pipe("cut -f1,55,56,110 /dcl01/leek/data/sraintrons/all_SRA_metadata.tsv")
depth = read.delim(cc, 	as.is=TRUE)
depth$totalReads = depth$sra_spots
pairIndex=which(depth$sra_layout == "PAIRED")
depth$totalReads[pairIndex] = 2*depth$sra_spots[pairIndex]

pheno$totalReads = depth$totalReads[match(pheno$sraRun, 
	depth$run_accession)]
pheno$totalReads10M = pheno$totalReads/10e6		

# i = Sys.getenv("SGE_TASK_ID") ## array job

for(i in seq(along=chrIndexes)) {
	### run script ##
	cat(".")
	ii = chrIndexes[[i]]
	calls = paste0("sed -n '",rowIndex$rowStart[i],",",
		rowIndex$rowEnd[i], "p; ", 
		rowIndex$rowEnd[i]+1, "q' ", fn, " | cut -f", 7:8)
		
	id = read.delim(pipe(calls[1]), header=FALSE,as.is=TRUE)$V1
	idList = strsplit(id, ",")

	cover =  read.delim(pipe(calls[2]), header=FALSE,as.is=TRUE)$V1
	coverList =  strsplit(cover, ",")

	dat = data.frame(row = rep(seq(along=idList),elementLengths(idList)),
		col = as.integer(unlist(idList)), 
		cover = as.integer(unlist(coverList)))
	dat$fullRow  = dat$row + rowIndex$skip[i]
	dat$colIndex= match(dat$col, pheno$colIndex)
	dat$depth = pheno$totalReads10M[
		match(dat$colIndex, pheno$colIndex)]
	dat$covAdj = dat$cover/dat$depth

	sm = sparseMatrix(i=dat$row, j = dat$colIndex, x = dat$cover,
		dims = c(length(id), nrow(pheno)),
		dimnames = list(names(jMap)[ii], pheno$sraRun))
	smadj = sparseMatrix(i=dat$row, j = dat$colIndex, 
		x = dat$covAdj,	dims = c(length(id), nrow(pheno)),
		dimnames = list(names(jMap)[ii], pheno$sraRun))
	assign(names(chrIndexes)[i], sm)
	assign(paste0(names(chrIndexes)[i],"_map"), jMap[ii])
	assign(paste0(names(chrIndexes)[i],"_adj"), smadj)

	save(list = c(names(chrIndexes)[i], 
			paste0(names(chrIndexes)[i],"_map"),
			paste0(names(chrIndexes)[i],"_adj")),
		file=paste0("byChr/", names(chrIndexes)[i],
			"_junctions.rda"))

}