# mkdir -p logs
# Rscript complete_chr_map.R > logs/complete_chr_map_log.txt 2>&1
library('devtools')

## Get chr info
chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('ucsc', 'length'))

## Read ucsc to ensembl mapping info
download.file('https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_UCSC2ensembl.txt', 'GRCh38_UCSC2ensembl.txt')
chrMapEns <- read.table('GRCh38_UCSC2ensembl.txt', header = FALSE, col.names = c('ucsc', 'ensembl'), colClasses = 'character', sep = '\t')

ens_195 <- match(chrInfo$ucsc, chrMapEns$ucsc)

## chrEBV is missing
chrInfo[is.na(ens_195), ]

## Add ensembl info
ens_194 <- ens_195[!is.na(ens_195)]
chrInfo$ensembl[!is.na(ens_195)] <- chrMapEns$ensembl[ens_194]

## Read ucsc to gencode mapping info
download.file('https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_ensembl2gencode.txt', 'GRCh38_ensembl2gencode.txt')
chrMapGen <- read.table('GRCh38_UCSC2ensembl.txt', header = FALSE, col.names = c('ucsc', 'gencode'), colClasses = 'character', sep = '\t')

gen_195 <- match(chrInfo$ucsc, chrMapGen$ucsc)

## Again, chrEBV is missing
chrInfo[is.na(gen_195), ]

gen_194 <- gen_195[!is.na(gen_195)]
chrInfo$gencode[!is.na(gen_195)] <- chrMapGen$gencode[gen_194]

write.table(chrInfo, file = 'hg38.ucsc.sizes.ensembl.gencode', quote = FALSE, sep = '\t', row.names = FALSE)

message(paste(Sys.time(), 'md5sum for GRCh38_UCSC2ensembl.txt'))
system('md5sum GRCh38_UCSC2ensembl.txt')
message(paste(Sys.time(), 'md5sum for GRCh38_ensembl2gencode.txt'))
system('md5sum GRCh38_ensembl2gencode.txt')

file.remove('GRCh38_UCSC2ensembl.txt')
file.remove('GRCh38_ensembl2gencode.txt')

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
devtools::session_info()
