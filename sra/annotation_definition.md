When we refer to "annotation" in the paper "Human splicing diversity across the Sequence Read Archive," we include most of the latest gene annotation tracks from the UCSC Genome Browser for hg19 and hg38 except in GENCODE's case, where latest releases were downloaded from GENCODE's website directly. We exclude algorithmic genome sequence-based predictions and Ensembl, which for the primary assembly is contained in GENCODE. All annotations were researched and downloaded on January 24, 2016. All hg38 exon-exon junction coordinates were lifted over to hg19 if that was possible.

We downloaded the files below and compressed them into our own archive "http://verve.webfactional.com/misc/jan_24_2016_annotations.tar.gz" to enable reproducibility.

### hg38 annotations

[UCSC genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz)

GENCODE v24 -- `ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz`

[RefSeq genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz)

[CCDS genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ccdsGene.txt.gz)

[MGC genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/mgcGenes.txt.gz)

[lincRNAs](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/lincRNAsTranscripts.txt.gz)

[SIB genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/sibGene.txt.gz)

### hg38 exceptions

Augustus (ab initio)
Geneid (ab initio)
Genscan (ab initio)
SGP (ab initio)

### hg19 annotations

[UCSC genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz)

[RefSeq genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz)

[AceView genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/acembly.txt.gz)

[CCDS genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ccdsGene.txt.gz)

GENCODE v19 -- `ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz`

[lincRNAs](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/lincRNAsTranscripts.txt.gz)

[MGC genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/mgcGenes.txt.gz)

[SIB genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/sibGene.txt.gz)

[Vega genes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vegaGene.txt.gz)

### hg19 exceptions

Augustus (ab initio)
Ensembl (contained in GENCODE)
Geneid (ab initio)
TransMap (cross-species)
Genscan (ab initio)
N-SCAN (ab initio)
