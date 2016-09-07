library(readr)
library(Matrix)
library(pbapply)
library(GenomicRanges)

getTiles <- function(counts, step=10000){
    p <- nrow(counts)
    tiles <- seq(1,p,step)
    if (!p %in% tiles){
        tiles <- c(tiles,p)
    }
    start <- tiles[-length(tiles)]
    end   <- tiles[-1]-1L
    tiles <- cbind(start=start, end=end)
    tiles
}

# Let's load the annotated junctions:
ann <- read.table("extdata/annotated_junctions.tsv")
ann.names <- paste0(ann$V1, ":", ann$V2, "-", ann$V3, "(", ann$V4, ")")


# Discarding samples with low coverage:
load("/dcl01/leek/data/sraintrons/filtered_sra_data_N1000_normCounts.rda")
good <- pd$totalReadsPred > 1e5
adjCounts <- adjCounts[,good]


# Annotated junctions:
ann.indices  <-  which(rownames(adjCounts) %in% ann.names) #248711 #new 248715
length(ann.indices)
adjCountsAnn <- adjCounts[ann.indices,]
totalJunctionCovAnn <- colSums(adjCountsAnn)
save(totalJunctionCovAnn, file="extdata/coverage_ann_info.rda")
tiles.ann    <- getTiles(adjCountsAnn)
save(tiles.ann, file="extdata/data/tiles_ann.rda")
for (i in 1:nrow(tiles.ann)){
    intM <- adjCountsAnn[tiles.ann[i,1]:tiles.ann[i,2],]
    intM <- intM - rowMeans(intM, na.rm=TRUE)
    save(intM, file=paste0("extdata/data/intM_ann_", i, ".rda"))
    print(i)
}
rm(adjCountsAnn)
gc()


# Unannotated junctions:
una.indices <-  which(!rownames(adjCounts) %in% ann.names)
length(una.indices) #56861
adjCountsUna <- adjCounts[una.indices,]
totalJunctionCovUna <- colSums(adjCountsUna)
save(totalJunctionCovUna, file="extdata/coverage_una_info.rda")
tiles.una <- getTiles(adjCountsUna)
save(tiles.una, file="extdata/data/tiles_una.rda")
for (i in 1:nrow(tiles.una)){
    intM <- adjCountsUna[tiles.una[i,1]:tiles.una[i,2],]
    intM <- intM - rowMeans(intM, na.rm=TRUE)
    save(intM, file=paste0("extdata/data/intM_una_", i, ".rda"))
    print(i)
}
rm(adjCountsUna)
gc()


# # All junctions:
# tiles <- getTiles(adjCounts)
# save(tiles, file="extdata/data/tiles_all.rda")
# totalJunctionCov <- colSums(adjCounts)
# totalNumberOfJunctions <- colSums(adjCounts!=0)
# njunctions <- nrow(adjCounts)
# save(totalJunctionCov, totalNumberOfJunctions,njunctions, file="extdata/coverage_info.rda")
# for (i in 1:nrow(tiles)){
#     intM <- adjCounts[tiles[i,1]:tiles[i,2],]
#     intM <- intM - rowMeans(intM, na.rm=TRUE)
#     save(intM, file=paste0("extdata/data/intM_all_", i, ".rda"))
#     print(i)
# }


