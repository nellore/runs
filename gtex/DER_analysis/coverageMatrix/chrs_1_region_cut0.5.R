chrs <- c("chrUn_KI270749v1"     ,   "chrUn_KI270435v1"
, "chr22_KI270732v1_random" ,"chrUn_KI270750v1"
 , "chr17_KI270730v1_random" ,"chr14_KI270724v1_random"
 , "chrUn_KI270519v1"        ,"chr14_KI270723v1_random"
 , "chrUn_KI270583v1"       , "chr9_KI270717v1_random"
, "chrUn_KI270589v1"        ,"chr5_GL000208v1_random"
, "chr14_GL000225v1_random" ,"chrM"
, "chr22_KI270735v1_random")

for(chr in chrs) {
    system(paste0('rm logs/*', chr, '*'))
    system(paste0('qsub .gtex-railMat-', chr, '.sh'))
}