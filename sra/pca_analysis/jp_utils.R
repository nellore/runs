
darken <- function(color, factor=1.4){
   col <- col2rgb(color)
   col <- col/factor
   col <- rgb(t(col), maxColorValue=255)
   col
}

getVar <- function(svd){
	d <- svd$d^2
	round(d/sum(d)*100,2)
}



getColorsBuci <- function(){
	col1  <- c(0.368417, 0.506779, 0.709798)
	col2  <- c(0.880722, 0.611041, 0.142051)
	col3  <- c(0.560181, 0.691569, 0.194885)
	col4  <- c(0.922526, 0.385626, 0.209179)
	col5  <- c(0.528488, 0.470624, 0.701351)
	col6  <- c(0.772079, 0.431554, 0.102387)
	col7  <- c(0.363898, 0.618501, 0.782349)
	col8  <- c(1, 0.75, 0)
	col9  <- c(0.647624, 0.37816, 0.614037)
	col10 <- c(0.571589, 0.586483, 0)
	col11 <- c(0.915, 0.3325, 0.2125)
	col12 <- c(0.40082222609352647, 0.5220066643438841, 0.85)
	col13 <- c(0.9728288904374106, 0.621644452187053, 0.07336199581899142)
	col14 <- c(0.736782672705901, 0.358, 0.5030266573755369)
	col15 <- c(0.28026441037696703, 0.715, 0.4292089322474965)
	colors <- rbind(col1,col2,col3,col4,col5,col6,col6,col8,col9,col10,col11,col12,col13,col14,col15)
	colors <- apply(colors,1,function(x){
		rgb(x[1], x[2], x[3])
	})
	colors
}


getMetadata <- function(indices=(1:21505-1)){
#### Now working with the metadata
	code <- read.table("metadata/index_to_SRA_accession.tsv", head=FALSE)
	#indices <- 1:21505-1
	code <- code[match(indices, code$V1),]
	thecall = pipe("cut -f1-88,90- /dcl01/leek/data/sraintrons/all_SRA_metadata.tsv")
	pd = read.delim(thecall, as.is=TRUE)
	rownames(pd) = pd$run_accession
	absCall = pipe("cut -f89 /dcl01/leek/data/sraintrons/all_SRA_metadata.tsv")
	ab = scan(absCall, sep="\n", what="character")
	pd$sra_study_abstract = ab[-1]
	metadata <- pd
	cc = pipe("cut -f1,55,56,110 /dcl01/leek/data/sraintrons/all_SRA_metadata.tsv")
	depth = read.delim(cc,     as.is=TRUE)
	depth$totalReads = depth$sra_spots
	pairIndex=which(depth$sra_layout == "PAIRED")
	depth$totalReads[pairIndex] = 2*depth$sra_spots[pairIndex]
	totalReads <- depth$totalReads
	totalReads10M <- totalReads/10e6 
	metadata$totalReads[match(depth$run_accession, metadata$run_accession)] <- totalReads 
	metadata$totalReads10M[match(depth$run_accession, metadata$run_accession)] <- totalReads10M 
	file <-  pipe("cut -d ',' -f1-7 metadata/sra-all-fields-2015-9-13.txt")
	tissue_info <-  read.delim(file, sep=",")
	indices <- match(metadata$run_accession, tissue_info$run)
	metadata$shark_tissue <- tissue_info$tissue[indices]
	metadata <- metadata[match(code$V5, metadata$run_accession),]
	metadata
}




getTissue <- function(metadata){
	metadata$shark
}




	