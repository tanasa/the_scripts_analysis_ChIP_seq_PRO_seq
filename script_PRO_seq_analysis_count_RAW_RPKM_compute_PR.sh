#!/bin/bash

################################################################################################

module load bowtie/2.3.4.2
module load homer/4.10

module load samtools/1.9
module load bedtools/2.26.0

module load fastx_toolkit/0.0.14
module load fastqc/0.11.7
module load trim_galore/0.4.5

################################################################################################

analyzeRNA.pl rna hg38 -count genes -strand + -pc 1 -noadj -start 200 -d \
PROseq.minusE2.bed.srm.dir/ \
PROseq.plusE2.bed.srm.dir/ \
> rna.hg38.PROseq.expression.RAW.start200.strand+.txt

analyzeRNA.pl rna hg38 -count genes -strand + -pc 1 -rpkm -start 200 -d \
PROseq.minusE2.bed.srm.dir/ \
PROseq.plusE2.bed.srm.dir/ \
> rna.hg38.PROseq.expression.RPKM.start200.strand+.txt

analyzeRepeats.pl rna hg38 -count pausing -strand + \
-pPromoterStart -50 -pPromoterEnd 200 -pBodyStart 200 -pBodyEnd 5000 -d \
PROseq.minusE2.bed.srm.dir/ \
PROseq.plusE2.bed.srm.dir/ \
> rna.hg38.PROseq.PAUSING.50_200.200_5000.strand+.txt

################################################################################################
