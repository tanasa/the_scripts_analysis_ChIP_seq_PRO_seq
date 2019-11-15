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
################################################################################################

annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
ChIRC_7SK_gRNA3_minusE2.Fan154_S21.Fan_154_S21_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
ChIRC_7SK_gRNA3_minusE2_FanC213.Fan_C213_S10_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
ChIRC_7SK_gRNA3_plusE2.Fan155_S22.Fan_155_S22_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.ChIRC_7SK_minusE2_plusE2.txt


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
ChIRC_7SK3_shctrl_minusE2_FanC117.Fan_C117_S11_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
ChIRC_7SK3_shHP1a1b_minusE2_FanC120.Fan_C120_S14_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
ChIRC_7SK3_shKAP1_minusE2_FanC119.Fan_C119_S13_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.ChIRC_7SK_schCTRL_shHP1_shKAP1.txt


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
G9a_minusE2_Fan83.Fan_83_S50_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
G9a_plusE2_Fan84.Fan_84_S51_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
G9a_minusE2_Fan91.Fan_91_S13_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
G9a_plusE2_Fan92.Fan_92_S14_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.G9a.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
H3K122ac_minusE2_Fan61.FAN_61_S39_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K122ac_plusE2_Fan62.FAN_62_S40_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.H3K122ac.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
H3K27me3_ctrlaso_minusE2_Fan25.Fan_25_S8_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K27me3_ctrlaso_plusE2_Fan26.Fan_26_S9_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K27me3_minusE2_Fan99.Fan_99_S7_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K27me3_plusE2_Fan100.Fan_100_S8_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.H3K27me3.txt       


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
H3K36me3_minusE2_Fan133.Fan_133_S51_L006_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K36me3_plusE2_Fan134.Fan_134_S52_L006_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K36me3_minusE2_Fan57.FAN_57_S87_L008_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K36me3_plusE2_Fan58.FAN_58_S88_L008_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.H3K36me3.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
H3K4me3_minusE2.Fan152_S19.Fan_152_S19_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K4me3_plusE2.Fan153_S20.Fan_153_S20_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K4me3_minusE2_Fan97.Fan_97_S5_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K4me3_plusE2_Fan98.Fan_98_S6_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.H3K4me3.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
H3K56ac_minusE2_Fan59.FAN_59_S37_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K56ac_plusE2_Fan60.FAN_60_S38_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K56ac_minusE2_Fan87.Fan_87_S54_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K56ac_plusE2_Fan88.Fan_88_S55_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.H3K56ac.txt 

annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
H3K9K14ac_minusE2_Fan135.Fan_135_S53_L006_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K9K14ac_plusE2_Fan136.Fan_136_S54_L006_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K9K14ac_minusE2_Fan85.Fan_85_S52_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K9K14ac_plusE2_Fan86.Fan_86_S53_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.H3K9K14ac.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
H3K9me3_minusE2_rep2_Fan104.Fan_ChIP_104_S4_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K9me3_plusE2_rep2_Fan103.Fan_ChIP_103_S3_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.H3K9me3.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
H3K9me3_shctrl_minusE2_FanC107.Fan_C107_S1_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K9me3_shKDM4B4C_minusE2_FanC109.Fan_C109_S3_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K9me3_shctrl_minusE2_FanC111.Fan_C111_S5_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
H3K9me3_shKDM4B4C_minusE2_FanC113.Fan_C113_S7_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.H3K9me3.schCTRL.shKDM4B.shKDM4C.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
HEXIM1_minusE2_FanC123.Fan_C123_S17_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
HEXIM1_plusE2_FanC124.Fan_C124_S18_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.HEXIM1.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
HP1_minusE2_Fan94.Fan_94_S2_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
HP1_plusE2_Fan93.Fan_93_S1_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.HP1.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
KAP1_minusE2_C140.Fan_C_140_S41_L007_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
KAP1_plusE2_C141.Fan_C_141_S42_L007_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
KAP1_minusE2_Fan137.Fan_137_S30_L007_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
KAP1_plusE2_Fan138.Fan_138_S31_L007_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.KAP1.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
KDM4B_minusE2_Fan89.Fan_89_S11_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
KDM4B_plusE2_Fan90.Fan_90_S12_L001_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
KDM4B_minusE2_FanC103.Fan_C103_S21_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
KDM4B_plusE2_FanC104.Fan_C104_S22_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.KDM4B.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
KDM4C_minusE2_Fan69.Fan_69_S47_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
KDM4C_plusE2_Fan70.Fan_70_S57_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
KDM4C_minusE2_FanC101.Fan_C101_S19_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
KDM4C_plusE2_FanC102.Fan_C102_S20_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.KDM4C.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
NELFA_minusE2.Fan150_S17.Fan_150_S17_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
NELFA_plusE2.Fan151_S18.Fan_151_S18_L005_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir / \
> tss.hg38.counts.NORM1e7.size.1kb.NELFA.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
NELFA_shctrl_minusE2_Fan114.Fan_ChIP_114_S39_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
NELFA_shctrl_plusE2_Fan113.Fan_ChIP_113_S38_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
NELFA_shHP1ab_minusE2_Fan116.Fan_ChIP_116_S41_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.NELFA.shCTRL.shHP1.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
POL2_ctrlaso_minusE2_Fan11.Fan_11_S15_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_plusE2_Fan12.Fan_12_S16_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_minusE2_Fan53.FAN53_S25_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
POL2_ctrlaso_plusE2_Fan54.FAN54_S26_L002_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.POL2.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
Pol2_shctrl_minusE2_Fan107.Fan_ChIP_107_S7_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
Pol2_shctrl_plusE2_Fan108.Fan_ChIP_108_S8_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
Pol2_shKDM4B4C_minusE2_Fan105.Fan_ChIP_105_S5_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
Pol2_shKDM4B4C_plusE2_Fan106.Fan_ChIP_106_S6_L003_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.POL2.shCTRL.shKDM4B.shKDM4C.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
Pol2S2P_shCtrl_minusE2_C135.Fan_C135_S53_L008_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
Pol2S2P_shCtrl_plusE2_C134.Fan_C134_S52_L008_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
Pol2S2P_shKAP1_minusE2_C136.Fan_C136_S54_L008_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
Pol2S2P_shKAP1_plusE2_C137.Fan_C137_S55_L008_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.POL2S2P.txt 


annotatePeaks.pl tss hg38 -size 1000 -norm 1e7 -strand both -d \
SUV39h1_minusE2_dsgf_Fan67.Fan_67_S45_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
SUV39h1_plusE2_dsgf_Fan68.Fan_68_S46_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
SUV39h1_minusE2_Fan66.Fan_66_S44_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
SUV39h1_plusE2_Fan65.Fan_65_S43_L004_R1_001.fastq.gz.after.trimmomatic.fastq.gz.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.SUV39h1.txt 


################################################################################################
################################################################################################
