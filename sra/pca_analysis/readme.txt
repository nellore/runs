Here are the different steps for the SVD:

Two input files are needed: 
a) annotated_junctions.tsv 
b) filtered_sra_data_N1000_normCounts.rda

1. Run splitJunctionsByBlock.R
2. Run createAta_ann.R and createAta_una.R through their shell scripts
3. Run createSVD.R
4. To create PCs for figure, run prepare_pca_data.R

# On the Hopkins cluster, the input files are at:
a) /dcl01/leek/data/gtex_work/runs/sra/annotated_junctions.tsv.gz
b) /dcl01/leek/data/sraintrons/filtered_sra_data_N1000_normCounts.rda