#pca.matrix.all <- read.table("/extdata/sra_junctions_pcmatrix_cutoff1000_cleaned.tsv", head=TRUE)
pca.matrix.ann <- read.table("extdata/sra_junctions_pcmatrix_cutoff1000_annotated.tsv", head=TRUE)
pca.matrix.una <- read.table("extdata/sra_junctions_pcmatrix_cutoff1000_unannotated.tsv", head=TRUE)
pca.matrices <- list(all=pca.matrix.all, ann=pca.matrix.ann, una=pca.matrix.una)
pca.matrices <- list(ann=pca.matrix.ann, una=pca.matrix.una)
pcs.list <- lapply(pca.matrices, function(pca.matrix){
	pcs <- lapply(1:10, function(i){ pca.matrix[,i]})
	pcs
})

source("jp_utils.R")
metadata <- getMetadata()
metadata <- metadata[match(rownames(pca.matrices[[1]]), metadata$run_accession),]
sra_tissue <- as.character(metadata$sra_tissue)
sra_source <- as.character(metadata$sra_source)
sra_cell   <- as.character(metadata$sra_cell_type)
sra_line   <- as.character(metadata$sra_cell_line)
sra_study  <- as.character(metadata$sra_study_accession)


# Let's make the sign of the pcs consistent across sets of junctions:
good <- is.finite(metadata$totalReads)
ss   <- sign(cor(metadata$totalReads[good], pcs.list[[1]][[1]][good]))
if (ss==-1){
	pcs.list[[1]][[1]] <- -pcs.list[[1]][[1]]
}
for (i in 1:10){
	sign <- sign(cor(pcs.list[[1]][[i]], pcs.list[[2]][[i]]))
	if (sign==-1){
		pcs.list[[2]][[i]] <- - pcs.list[[2]][[i]]
	}
	#sign <- sign(cor(pcs.list[[1]][[i]], pcs.list[[3]][[i]]))
	#if (sign==-1){
	#		pcs.list[[3]][[i]] <- - pcs.list[[3]][[i]]
	#}
}


# pcs for unannotated junctions:
pcs <- do.call(cbind,pcs.list$una)
pcs <- as.data.frame(pcs[,1:5])
names(pcs) <- paste0("PC", 1:5)

# Brain indices:
indices1 <- which(grepl(tolower("BRAIN|PITUITARY|CORTEX|GLIOMA|GLIOBLASTOMA|NEURON|NEURAL|NEURONS"),tolower(sra_tissue)))
indices2 <- which(grepl(tolower("BRAIN|PITUITARY|CORTEX|GLIOBLASTOMA|GLIOMA|NEURON|NEURAL|NEURONS"),tolower(sra_cell)))
indices3 <- which(grepl(tolower("BRAIN|PITUITARY|CORTEX|GLIOBLASTOMA|HBRR|GLIOMA|NEURON|NEURAL|NEURONS"),tolower(sra_source)))
indices.abrf <- which(sra_study=="SRP026126")
indices.seqc <- which(sra_study=="SRP025982")
indices.brain <- unique(c(indices1,indices2,indices3))
indices.brain <- setdiff(indices.brain, indices.seqc)

# Blood indices:
indices1 <- which(grepl(tolower("BLOOD|MONOCYTE|PERIPHERAL"),tolower(sra_tissue)))
indices2 <- which(grepl(tolower("BLOOD|MONOCYTE|PERIPHERAL"),tolower(sra_cell)))
indices3 <- which(grepl(tolower("BLOOD|MONOCYTE|PERIPHERAL"),tolower(sra_source)))
indices.blood <- unique(c(indices1,indices2,indices3))

# LCL indicesL
indices1 <- which(sra_study=="SRP001563")
indices2 <- which(sra_study=="SRP001540")
indices3 <- which(sra_study=="SRP030041")
indices.geuvadis <- which(sra_study=="ERP001942")
indices.lcl <- unique(c(indices1,indices2,indices3, indices.geuvadis))

# Adding to the data frame:
pcs$blood <- pcs$lcl <- pcs$brain <- pcs$seqc_abrf <- pcs$geuvadis <- rep(0, nrow(pcs))
pcs$blood[indices.blood] <- 1
pcs$brain[indices.brain] <- 1
pcs$lcl[indices.lcl] <- 1
pcs$seqc_abrf[indices.seqc] <- 1
pcs$seqc_abrf[indices.abrf] <- 2
pcs$geuvadis[indices.geuvadis] <- 1

# Adding the ratios for SEQC:
# Note that we are inversing ratios so that HBRR is first
pcs$mixture <- rep(NA, nrow(pcs))
pcs$mixture[indices.seqc[sra_source[indices.seqc]=="UNIVERSAL HUMAN REFERENCE RNA (UHRR) FROM STRATAGENE"]] <- "0:1"
pcs$mixture[indices.seqc[sra_source[indices.seqc]=="HUMAN BRAIN REFERENCE RNA (HBRR) FROM AMBION"]] <- "1:0"
pcs$mixture[indices.seqc[sra_source[indices.seqc]=="MIX OF SAMPLE A (UHRR) AND SAMPLE B (HBRR) (1:3)"]] <- "3:1"
pcs$mixture[indices.seqc[sra_source[indices.seqc]=="MIX OF SAMPLE A (UHRR) AND SAMPLE B (HBRR) (3:1)"]] <- "1:3"
pcs$mixture[indices.seqc[sra_source[indices.seqc]=="MIX OF A AND B (1:3)"]] <- "3:1"
pcs$mixture[indices.seqc[sra_source[indices.seqc]=="MIX OF A AND B (3:1)"]] <- "1:3"


# Adding the ratios for ABRF:
# A: UHRR
# B: HBRR
# C: UHRR:HBRR = 3:1
# D: UHRR:HBRR = 1:3
labels <- strsplit(sra_source[indices.abrf], "-")
labels <- unlist(lapply(labels, function(x){
	n <- length(x)
	x[[n]]
}))
labels <- gsub("AH|AS|AR", "A",labels)
labels <- gsub("BR", "B",labels)
pcs$mixture[indices.abrf[labels=="A"]] <- "0:1"
pcs$mixture[indices.abrf[labels=="B"]] <- "1:0"
pcs$mixture[indices.abrf[labels=="C"]] <- "1:3"
pcs$mixture[indices.abrf[labels=="D"]] <- "3:1"
rownames(pcs) <- rownames(pca.matrices[[1]])

write.table(pcs, row.names=TRUE, col.names=TRUE, quote=FALSE, file="../pcs_unannotated_with_pd.tsv")



# Calculate correlations:
load("../extdata/coverage_una_info.rda")
load("../extdata/coverage_info.rda")
#load("../extdata/coverage_raw_info.rda")


#good <- match(rownames(pcs),names(totalJunctionCovRaw))
#cor(log2(totalJunctionCovRaw[good]+1), pcs[,1]) 


# Correlation with log2(mapped reads)
good <- is.finite(metadata$totalReads) 
cor(log2(metadata$totalReads[good]+1), pcs[good,1]) #0.277


# Correlation with read length
pd <- metadata
pd$readLength = pd$sra_bases / pd$sra_spots
pd$readLength[which(pd$sra_layout == "PAIRED")] = pd$readLength[which(pd$sra_layout == "PAIRED")]/2
pd$readLength = round(pd$readLength)
good <- is.finite(pd$readLength) 
cor(log2(pd$readLength[good]+1), pcs[good,1]) #0.639

cor(totalJunctionCovUna, pcs[,1]) #0.9768495

#cor(totalNumberOfJunctions, pcs[,1])
#cor(totalJunctionCov, pcs[,1])







