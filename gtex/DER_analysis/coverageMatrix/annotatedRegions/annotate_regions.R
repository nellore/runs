library('GenomicRanges')
library('derfinder')
library('devtools')
library('GenomeInfoDb')
library('derfinderPlot')
library('limma')

## Load chr matching info
chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/genomicState/hg38.ucsc.sizes.ensembl.gencode', header = TRUE, colClasses = c('character', 'numeric', 'character', 'character'))

## Load regions
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/regions-cut0.5.Rdata')

## Format seqnames in Ensembl format, drop any chr that doesn't match (it's only chrEBV)
regs_ensembl <- regions

ensembl <- chrInfo$ensembl[match(seqlevels(regs_ensembl), chrInfo$ucsc)]
regs_ensembl <- regs_ensembl[!seqnames(regs_ensembl) %in% seqlevels(regs_ensembl)[is.na(ensembl)]]

regs_ensembl <- dropSeqlevels(regs_ensembl, seqlevels(regs_ensembl)[is.na(ensembl)])
seqlevels(regs_ensembl) <- ensembl[!is.na(ensembl)]

## Annotate regions by ensembl
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/genomicState/genomicState.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5.Rdata')
annotated_ensembl <- annotateRegions(regions = regs_ensembl, genomicState = genomicState.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5$fullGenome, chrsStyle = NULL)
annotated_ensembl_8 <- annotateRegions(regions = regs_ensembl, genomicState = genomicState.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5$fullGenome, chrsStyle = NULL, minoverlap = 8L)
annotated_ensembl_20 <- annotateRegions(regions = regs_ensembl, genomicState = genomicState.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5$fullGenome, chrsStyle = NULL, minoverlap = 20L)
save(annotated_ensembl, file = 'annotated_ensembl.Rdata')
save(annotated_ensembl_8, file = 'annotated_ensembl_8.Rdata')
save(annotated_ensembl_20, file = 'annotated_ensembl_20.Rdata')
save(regs_ensembl, file = 'regs_ensembl.Rdata')

## Annotate regions by UCSC -- dropping chrEBV because it's not on the genomic state object
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/genomicState/genomicState.Hsapiens.UCSC.hg38.knownGene.Rdata')
regs_ucsc <- dropSeqlevels(regions, 'chrEBV')
message(paste(Sys.time(), 'there are', length(regions) - length(regs_ucsc), 'regions in chrEBV that are excluded'))
annotated_ucsc <- annotateRegions(regions = regs_ucsc, genomicState = genomicState.Hsapiens.UCSC.hg38.knownGene$fullGenome, chrsStyle = NULL)
annotated_ucsc_8 <- annotateRegions(regions = regs_ucsc, genomicState = genomicState.Hsapiens.UCSC.hg38.knownGene$fullGenome, chrsStyle = NULL, minoverlap = 8L)
annotated_ucsc_20 <- annotateRegions(regions = regs_ucsc, genomicState = genomicState.Hsapiens.UCSC.hg38.knownGene$fullGenome, chrsStyle = NULL, minoverlap = 20L)
save(annotated_ucsc, file = 'annotated_ucsc.Rdata')
save(annotated_ucsc_8, file = 'annotated_ucsc_8.Rdata')
save(annotated_ucsc_20, file = 'annotated_ucsc_20.Rdata')
save(regs_ucsc, file = 'regs_ucsc.Rdata')


vennInfoUc <- vennInfoEns <- list(
    'min1bp' = list('1bp' = NA, '8bp' = NA, '20bp' = NA),
    'min8bp' = list('1bp' = NA, '8bp' = NA, '20bp' = NA),
    'min20bp' = list('1bp' = NA, '8bp' = NA, '20bp' = NA)
)
## Make venn plots -- min 1 bp overlap
pdf('venn_regions_min1bp_overlap.pdf')
vennInfoEns[['min1bp']][['1bp']] <- vennRegions(annotated_ensembl, main = 'Regions by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
vennInfoUc[['min1bp']][['1bp']] <- vennRegions(annotated_ucsc, main = 'Regions by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 7 bp
vennInfoEns[['min1bp']][['8bp']] <- vennRegions(annotated_ensembl, subsetIndex = width(regs_ensembl) > 7, main = 'Regions > 7 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
vennInfoUc[['min1bp']][['8bp']] <- vennRegions(annotated_ucsc, subsetIndex = width(regs_ucsc) > 7, main = 'Regions > 7 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 19 bp
vennInfoEns[['min1bp']][['20bp']] <- vennRegions(annotated_ensembl, subsetIndex = width(regs_ensembl) > 19, main = 'Regions > 19 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
vennInfoUc[['min1bp']][['20bp']] <- vennRegions(annotated_ucsc, subsetIndex = width(regs_ucsc) > 19, main = 'Regions > 19 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')
dev.off()

## Make venn plots -- min 9 bp overlap
pdf('venn_regions_min8bp_overlap.pdf')
vennInfoEns[['min8bp']][['1bp']] <- vennRegions(annotated_ensembl_8, main = 'Regions by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
vennInfoUc[['min8bp']][['1bp']] <- vennRegions(annotated_ucsc_8, main = 'Regions by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 7 bp
vennInfoEns[['min8bp']][['8bp']] <- vennRegions(annotated_ensembl_8, subsetIndex = width(regs_ensembl) > 7, main = 'Regions > 7 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
vennInfoUc[['min8bp']][['8bp']] <- vennRegions(annotated_ucsc_8, subsetIndex = width(regs_ucsc) > 7, main = 'Regions > 7 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 19 bp
vennInfoEns[['min8bp']][['20bp']] <- vennRegions(annotated_ensembl_8, subsetIndex = width(regs_ensembl) > 19, main = 'Regions > 19 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
vennInfoUc[['min8bp']][['20bp']] <- vennRegions(annotated_ucsc_8, subsetIndex = width(regs_ucsc) > 19, main = 'Regions > 19 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')
dev.off()


## Make venn plots -- min 20 bp overlap
pdf('venn_regions_min20bp_overlap.pdf')
vennInfoEns[['min20bp']][['1bp']] <- vennRegions(annotated_ensembl_20, main = 'Regions by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
vennInfoUc[['min20bp']][['1bp']] <- vennRegions(annotated_ucsc_20, main = 'Regions by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 7 bp
vennInfoEns[['min20bp']][['8bp']] <- vennRegions(annotated_ensembl_20, subsetIndex = width(regs_ensembl) > 7, main = 'Regions > 7 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
vennInfoUc[['min20bp']][['8bp']] <- vennRegions(annotated_ucsc_20, subsetIndex = width(regs_ucsc) > 7, main = 'Regions > 7 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 19 bp
vennInfoEns[['min20bp']][['20bp']] <- vennRegions(annotated_ensembl_20, subsetIndex = width(regs_ensembl) > 19, main = 'Regions > 19 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
vennInfoUc[['min20bp']][['20bp']] <- vennRegions(annotated_ucsc_20, subsetIndex = width(regs_ucsc) > 19, main = 'Regions > 19 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')
dev.off()

vennInfo <- list('ensembl' = vennInfoEns, 'ucsc' = vennInfoUc)

vennInfo <- lapply(vennInfo, function(x) { 
    lapply(x, function(y) {
        lapply(y, function (z) {
            res <- cbind(z, 'Percent' = z[,'Counts'] * 100 / sum(z[, 'Counts']))
            return(res)
        })
    })
})

plotPercent <- function(x, ...) {
    info <- structure(cbind(x[, 1:3], Counts= round(x[, 'Percent'], 2)), class="VennCounts")
    vennDiagram(info, ...)
}

## Percent min overlap 1 bp
pdf('venn_regions_min1bp_overlap_percent.pdf')
plotPercent(vennInfo[['ensembl']][['min1bp']][['1bp']], main = 'Regions by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
plotPercent(vennInfo[['ucsc']][['min1bp']][['1bp']], main = 'Regions by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 7 bp
plotPercent(vennInfo[['ensembl']][['min1bp']][['8bp']], main = 'Regions > 7 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
plotPercent(vennInfo[['ucsc']][['min1bp']][['8bp']], main = 'Regions > 7 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 19 bp
plotPercent(vennInfo[['ensembl']][['min1bp']][['20bp']], main = 'Regions > 19 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
plotPercent(vennInfo[['ucsc']][['min1bp']][['20bp']], main = 'Regions > 19 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')
dev.off()

## Percent min overlap 8 bp
pdf('venn_regions_min8bp_overlap_percent.pdf')
plotPercent(vennInfo[['ensembl']][['min8bp']][['1bp']], main = 'Regions by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
plotPercent(vennInfo[['ucsc']][['min8bp']][['1bp']], main = 'Regions by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 7 bp
plotPercent(vennInfo[['ensembl']][['min8bp']][['8bp']], main = 'Regions > 7 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
plotPercent(vennInfo[['ucsc']][['min8bp']][['8bp']], main = 'Regions > 7 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 19 bp
plotPercent(vennInfo[['ensembl']][['min8bp']][['20bp']], main = 'Regions > 19 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
plotPercent(vennInfo[['ucsc']][['min8bp']][['20bp']], main = 'Regions > 19 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')
dev.off()

## Percent min overlap 20 bp
pdf('venn_regions_min20bp_overlap_percent.pdf')
plotPercent(vennInfo[['ensembl']][['min20bp']][['1bp']], main = 'Regions by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
plotPercent(vennInfo[['ucsc']][['min20bp']][['1bp']], main = 'Regions by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 7 bp
plotPercent(vennInfo[['ensembl']][['min20bp']][['8bp']], main = 'Regions > 7 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
plotPercent(vennInfo[['ucsc']][['min20bp']][['8bp']], main = 'Regions > 7 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')

## Greater than 19 bp
plotPercent(vennInfo[['ensembl']][['min20bp']][['20bp']], main = 'Regions > 19 bp by ENSEMBL GRCh38.p5 annotation (excluding chrEBV)', counts.col = 'blue')
plotPercent(vennInfo[['ucsc']][['min20bp']][['20bp']], main = 'Regions > 19 bp by UCSC knownGene hg38 annotation (excluding chrEBV)', counts.col = 'blue')
dev.off()

## Save information
save(vennInfo, file = 'vennInfo.Rdata')

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
