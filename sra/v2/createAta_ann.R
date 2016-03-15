i=as.numeric(commandArgs(TRUE)[1]) # Index of the tile
load("extdata/data/tiles_ann.rda")
n <- nrow(tiles.ann)
files <- paste0("extdata/data/intM_ann_",1:n, ".rda")
file  <- files[i]

cat("Loading data... \n")
load(file)

cat("Crossproduct...")
ata <- crossprod(intM)

cat("Saving data... \n")
outputfile <- gsub("intM", "ata", file)
save(ata, file=outputfile)
q("no")






