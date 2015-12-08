## This script reads the phenotype data

pheno <- read.table('/dcl01/leek/data/gtex_work/gtex_pheno_table.tsv', sep = '\t', header = TRUE, comment.char = '', na.strings = c('NA', 'N/A'), quote = '')
cols.96.97.98.99 <- c('RACE', 'TRAMP', 'DTHATPSY', 'DTHRFG', 'DTHCERT', 'DTHVNT', 'DTHWTNS', 'LBCMVTAB', 'LBEBVGAB', 'LBEBVMAB', 'LBHBCABM', 'LBHBCABT', 'LBHBSAB', 'LBHBSAG', 'LBHCV1NT', 'LBHBHCVAB', 'LBHIV1NT', 'LBHIVAB', 'LBHIVO', 'LBPRRVDRL', 'LBRPR', 'MHABNWBC', 'MHALS', 'MHALZDMT', 'MHALZHMR', 'MHARTHTS', 'MHASCITES', 'MHASTHMA', 'MHBCTINF', 'MHBLDDND', 'MHBLDOCNT', 'MHCANCER5', 'MHCANCERC', 'MHCANCERNM', 'MHCLLULTS', 'MHCLRD', 'MHCOCAINE5', 'MHCOPD', 'MHCOUGHU', 'MHCVD', 'MHDLYSIS', 'MHDMNTIA', 'MHDPRSSN', 'MHDTND72H', 'MHENCEPHA', 'MHEURO5', 'MHFLU', 'MHFNGINF', 'MHFVRU', 'MHGNRR12M', 'MHHEPBCT', 'MHHEPCCT', 'MHHEROIN', 'MHHGH', 'MHHIVCT', 'MHHIVNT', 'MHHMPHLIA', 'MHHMPHLIAB', 'MHHRTATT', 'MHHRTDIS', 'MHHRTDISB', 'MHHTN', 'MHINFLNE', 'MHIVDRG5', 'MHJAKOB', 'MHLAPTHU', 'MHLUPUS', 'MHLVRDIS', 'MHMENINA', 'MHMS', 'MHMSXWMA', 'MHMSXWMB', 'MHNEPH', 'MHNGHTSWT', 'MHNPHYS4W', 'MHNRTHEUR', 'MHOPNWND', 'MHOPPINF', 'MHORGNTP', 'MHOSTMYLTS', 'MHPLLABS', 'MHPNMIAB', 'MHPNMNIA', 'MHPRCNP', 'MHPRKNSN', 'MHPSBLDCLT', 'MHRA', 'MHRBSANML', 'MHREYES', 'MHRNLFLR', 'MHSARS', 'MHSCHZ', 'MHSCLRDRM', 'MHSDRGABS', 'MHSEPSIS', 'MHSKNSPT', 'MHSMLPXCT', 'MHSMLPXVC', 'MHSRC', 'MHSRCDSS', 'MHSRGHM', 'MHSTD', 'MHSTRDLT', 'MHSUBABSA', 'MHSUBABSB', 'MHSXMDA', 'MHSXMDB', 'MHSYPH12M', 'MHSZRSU', 'MHT1D', 'MHT2D', 'MHTBHX', 'MHTEMPU', 'MHTTOO12M', 'MHTTOONP', 'MHTXCEXP', 'MHUK8096', 'MHUREMIA', 'MHWKNSSU', 'MHWNVCT', 'MHWNVHX', 'MHWTLSUA', 'MHWTLSUB')
stopifnot(all(cols.96.97.98.99 %in% colnames(pheno)))
stopifnot(length(cols.96.97.98.99) == 109 + 14)

for(i in cols.96.97.98.99) pheno[pheno[, i] %in% 96:99, i] <- NA

## Label GENDER
pheno$GENDER <- factor(pheno$GENDER, levels = 1:2, labels = c('male', 'female'))
stopifnot(all(pheno$GENDER == pheno$Sex))

## Label RACE
pheno$RACE <- factor(pheno$RACE, levels = 1:4, labels = c('Asian', 'Black or African American', 'White', 'American Indian or Alaska Native'))

## Label ETHNCTY
pheno$ETHNCTY <- factor(pheno$ETHNCTY, levels = 0:1, labels = c('Not Hispanic or Latino', 'Hispanic or Latino'))

## Core body temperarture in celcius
pheno$TRCRTMP_Celcius <- pheno$TRCRTMP
pheno$TRCRTMP_Celcius[pheno$TRCRTMPU == 'F' & !is.na(pheno$TRCRTMP)] <- (pheno$TRCRTMP_Celcius[pheno$TRCRTMPU == 'F' & !is.na(pheno$TRCRTMP)] - 32) / 1.8

pheno$DTHHRDY <- factor(pheno$DTHHRDY, levels = 0:4, labels = c('Ventilator Case', 'Violent and fast death', 'Fast death of natural causes', 'Intermediate death', 'Slow Death'))

cols.neg.pos <- c('LBCMVTAB', 'LBEBVGAB', 'LBEBVMAB', 'LBHBCABM', 'LBHBCABT', 'LBHBSAB', 'LBHBSAG', 'LBHCV1NT', 'LBHBHCVAB', 'LBHIV1NT', 'LBHIVAB', 'LBHIVO', 'LBPRRVDRL', 'LBRPR')
stopifnot(length(cols.neg.pos) == 14)
for(i in cols.neg.pos) pheno[, i] <- factor(pheno[, i], levels = 0:1, labels = c('Negative', 'Positive'))

stopifnot(length(cols.96.97.98.99[ !cols.96.97.98.99 %in% c(cols.neg.pos, 'RACE', 'GENDER') ]) == 108)
for(i in cols.96.97.98.99[ !cols.96.97.98.99 %in% c(cols.neg.pos, 'RACE', 'GENDER') ]) pheno[, i] <- factor(pheno[, 1], levels = 0:1, labels = c('No', 'Yes'))

## Format ReleaseDate and LoadDate
for(i in c('ReleaseDate', 'LoadDate')) pheno[, i] <- as.Date(as.character(pheno[, i]), format = c('%b %d, %Y'))

## Format SMNABTCHD and SMGEBTCHD
for(i in c('SMNABTCHD', 'SMGEBTCHD')) pheno[, i] <- as.Date(as.character(pheno[, i]), format = c('%m/%d/%Y'))
    
## Format DTHPRNINT in minutes
pheno$DTHPRNINT_Minute <- as.integer(gsub(' hour.*', '', as.character(pheno$DTHPRNINT))) * 60 + as.integer(gsub('.*hour\\(s\\), | minute.*', '', as.character(pheno$DTHPRNINT)))

## Note that Rail RNA Batch number has to be a factor
pheno$RailRnaBatchNumber <- as.factor(pheno$RailRnaBatchNumber)


## Find which variables are all NA
number_NA <- sapply(pheno, function(x) { sum(is.na(x))})
names(number_NA) <- colnames(pheno)

## Number of variables with all NAs and their column names
all_NA <- colnames(pheno)[number_NA == nrow(pheno)]
length(all_NA)
all_NA

## Percent of missing observations on the remaining variables
summary(number_NA[!names(number_NA) %in% all_NA] / nrow(pheno) * 100)

## Variables with more than 10% of the observations missing (excluding those with 100%)
missing_10 <- colnames(pheno)[number_NA / nrow(pheno) > 0.1 & number_NA != nrow(pheno)]
summary(number_NA[missing_10] / nrow(pheno) * 100)
length(missing_10)

## Explore the variables missing more than 10% of the data (but not 100%)
summary(pheno[, missing_10])

## Explore data from variables missing between (0, 10]% of the data
sum(number_NA / nrow(pheno) <= 0.1)
sum(number_NA / nrow(pheno) <= 0.1 & number_NA / nrow(pheno) > 0)
summary(pheno[, number_NA / nrow(pheno) <= 0.1 & number_NA / nrow(pheno) > 0])

## Variables not missing anything
sum(number_NA / nrow(pheno) == 0)
summary(pheno[, number_NA / nrow(pheno) == 0])

save(pheno, file = 'pheno_complete.Rdata')
save(pheno[, number_NA / nrow(pheno) <= 0.1], file = 'pheno_missing_less_10.Rdata')


## Reproducibility info
Sys.time()
proc.time()
devtools::session_info()
