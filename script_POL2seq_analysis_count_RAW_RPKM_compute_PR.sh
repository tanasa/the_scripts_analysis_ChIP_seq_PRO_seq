#!/bin/bash

################################################################################################
################################################################################################

module load legacy/scg4

module load bowtie/2.3.4.2
module load homer/4.10

module load samtools/1.9
module load bedtools/2.26.0

module load fastx_toolkit/0.0.14
module load fastqc/0.11.7
module load trim_galore/0.4.5

################################################################################################

analyzeRNA.pl rna hg38 -count genes -strand both -pc 1 -noadj -start 200 -d \
POL2_ctrlaso_minusE2_Fan11.Fan_11_S15_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_plusE2_Fan12.Fan_12_S16_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_minusE2_Fan53.FAN53_S25_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_plusE2_Fan54.FAN54_S26_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> rna.hg38.POL2_seq.expression.RAW.start200.strand_both.txt

analyzeRNA.pl rna hg38 -count genes -strand both -pc 1 -rpkm -start 200 -d \
POL2_ctrlaso_minusE2_Fan11.Fan_11_S15_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_plusE2_Fan12.Fan_12_S16_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_minusE2_Fan53.FAN53_S25_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_plusE2_Fan54.FAN54_S26_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> rna.hg38.POL2_seq.expression.RPKM.start200.strand_both.txt

analyzeRepeats.pl rna hg38 -count pausing -strand both \
-pPromoterStart -50 -pPromoterEnd 200 -pBodyStart 200 -pBodyEnd 5000 -d \
POL2_ctrlaso_minusE2_Fan11.Fan_11_S15_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_plusE2_Fan12.Fan_12_S16_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_minusE2_Fan53.FAN53_S25_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_plusE2_Fan54.FAN54_S26_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> rna.hg38.POL2_seq.PAUSING.50_200.200_5000.strand_both.txt

################################################################################################
################################################################################################
