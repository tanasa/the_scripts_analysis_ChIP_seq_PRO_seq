#!/bin/bash

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

analyzeRNA.pl rna hg38 -count genes -strand + -pc 1 -noadj -start 200 -d \
PROseq.R32_shCTRL_rep2.RM032.RM-032_S4_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R30_shHP1ab_rep2.RM030.RM-030_S2_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R33_shCTRL_new.RM033.RM-033_S5_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R34_shHP1ab_new.RM034.RM-034_S6_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R37_shCTRL_minusE2.RM037.RM-037_S9_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R38_shCTRL_plusE2.RM038.RM-038_S10_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R39_shKAP1_minusE2.RM039.RM-039_S11_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R40_shKAP1_plusE2.RM040.RM-040_S12_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> rna.hg38.PROseq.expression.RAW.start200.strand+.txt

analyzeRNA.pl rna hg38 -count genes -strand + -pc 1 -rpkm -start 200 -d \
PROseq.R32_shCTRL_rep2.RM032.RM-032_S4_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R30_shHP1ab_rep2.RM030.RM-030_S2_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R33_shCTRL_new.RM033.RM-033_S5_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R34_shHP1ab_new.RM034.RM-034_S6_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R37_shCTRL_minusE2.RM037.RM-037_S9_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R38_shCTRL_plusE2.RM038.RM-038_S10_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R39_shKAP1_minusE2.RM039.RM-039_S11_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R40_shKAP1_plusE2.RM040.RM-040_S12_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> rna.hg38.PROseq.expression.RPKM.start200.strand+.txt

analyzeRepeats.pl rna hg38 -count pausing -strand + \
-pPromoterStart -50 -pPromoterEnd 200 -pBodyStart 200 -pBodyEnd 5000 -d \
PROseq.R32_shCTRL_rep2.RM032.RM-032_S4_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R30_shHP1ab_rep2.RM030.RM-030_S2_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R33_shCTRL_new.RM033.RM-033_S5_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R34_shHP1ab_new.RM034.RM-034_S6_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R37_shCTRL_minusE2.RM037.RM-037_S9_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R38_shCTRL_plusE2.RM038.RM-038_S10_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R39_shKAP1_minusE2.RM039.RM-039_S11_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
PROseq.R40_shKAP1_plusE2.RM040.RM-040_S12_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> rna.hg38.PROseq.PAUSING.50_200.200_5000.strand+.txt

################################################################################################
