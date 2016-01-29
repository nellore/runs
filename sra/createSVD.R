library(fsvd)

# For Filter n>=1000
ata.list <- list()
n <- 31
files <- paste0(paste0("extdata/data/ata_cleaned_",1:n, ".rda"))
load(files[1])
ata_full <- matrix(0, nrow(ata), nrow(ata))
# Loop to avoid memory complains:
for (i in 1:n){
	load(files[i])
	ata_full <- ata_full + ata
	print(i)
}
ata <- ata_full
save(ata, file="extdata/ata_cleaned.rda")


svd <- fsvd(ata, k=100)
pca.matrix <- svd$v
rownames(pca.matrix) <- colnames(ata)
colnames(pca.matrix) <- paste0("PC", 1:100)
write.table(pca.matrix, row.names=TRUE, col.names=TRUE, 
	quote=FALSE, file="extdata/sra_junctions_pcmatrix_cutoff1000_cleaned.tsv")



# For Filter n>=1000
# ANNOTATION:
ata.list <- list()
n <- 25
files <- paste0(paste0("extdata/data/ata_ann_",1:n, ".rda"))
load(files[1])
ata_full <- matrix(0, nrow(ata), nrow(ata))

for (i in 1:n){
	load(files[i])
	ata_full <- ata_full + ata	
	print(i)
}
ata <- ata_full
save(ata, file="extdata/ata_ann.rda")

svd <- fsvd(ata, k=100)
pca.matrix <- svd$v
rownames(pca.matrix) <- colnames(ata)
colnames(pca.matrix) <- paste0("PC", 1:100)
write.table(pca.matrix, row.names=TRUE, col.names=TRUE, 
	quote=FALSE, file="extdata/sra_junctions_pcmatrix_cutoff1000_annotated.tsv")



# For Filter n>=1000
# NO ANNOTATION:
ata.list <- list()
n <- 6
files <- paste0(paste0("extdata/data/ata_una_",1:n, ".rda"))
load(files[1])
ata_full <- matrix(0, nrow(ata), nrow(ata))

for (i in 1:n){
	load(files[i])
	ata_full <- ata_full + ata
	print(i)
}
ata <- ata_full
save(ata, file="extdata/ata_una.rda")

svd <- fsvd(ata, k=100)
pca.matrix <- svd$v
rownames(pca.matrix) <- colnames(ata)
colnames(pca.matrix) <- paste0("PC", 1:100)
write.table(pca.matrix, row.names=TRUE, col.names=TRUE, 
	quote=FALSE, file="extdata/sra_junctions_pcmatrix_cutoff1000_unannotated.tsv")




getVar <- function(svd){
	round(svd$d^2/sum(svd$d^2)*100,2)
}

#getVar(svd)[1:10]







