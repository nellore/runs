#!/usr/bin/env bash
# Runs Rail-RNA on translatome samples from project http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP031883
# Use Rail-RNA v0.2.3b and mm10 iGenome from ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
# Specify path to Bowtie 1/2 index basenames for mm10 iGenome here
MM10=/scratch2/langmead-fs1/mm10/Mus_musculus/UCSC/mm10/Sequence/
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna go local -m $DIR/translatome.manifest -x $MM10/BowtieIndex/genome,$MM10/Bowtie2Index/genome -d jx
