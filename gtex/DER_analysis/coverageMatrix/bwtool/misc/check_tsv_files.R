# qrsh -l mem_free=10G,h_vmem=11G -pe local 25

library('R.utils')
library('BiocParallel')

## Parallel environment
bp <- MulticoreParam(workers = 25, outfile = Sys.getenv('SGE_STDERR_PATH'))

## Load pheno data
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/pheno/pheno_missing_less_10.Rdata')

## Locate bwtool results
tsv <- dir('/dcl01/leek/data/gtex_work/mean_cov', pattern = 'tsv', full.names = TRUE)
names(tsv) <- gsub('.sum.tsv', '', dir('/dcl01/leek/data/gtex_work/mean_cov/', pattern = 'tsv'))

## Find which samples have tsv files and match them that way
i <- match(as.character(pheno$Run), names(tsv))
if(any(is.na(i))) message(paste(Sys.time(), 'there are', sum(is.na(i)), 'tsv files missing'))
## Ok, all files are there
stopifnot(length(tsv) == nrow(pheno))

## Sort tsv files to match pheno data
tsv2 <- tsv[i[!is.na(i)]]
j <- match(names(tsv2), as.character(pheno$Run))
identical(j, seq_len(length(j)))

sampleNames <- gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', pheno$BigWigPath[!is.na(i)])

## Check that all files have the same number of rows
system.time( tsv_lines <- bplapply(tsv[i], countLines, BPPARAM = bp) )
## Ok, they all have the same number of rows
table(unlist(tsv_lines))


suspect <- 'SRR2167642_SRS1036572_SRX1153642_male_brain.cerebellar.hemisphere'
which(gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', pheno$BigWigPath) == suspect)
suspect_regs <- c(447134, 447137, 447144)
suspect_tsv <- tsv2[4]
suspect_tsv

for(i in suspect_regs) system(paste('head -n', i, suspect_tsv, '| tail -n 1'))







