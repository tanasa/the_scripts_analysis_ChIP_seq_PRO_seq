#!/bin/bash

################################################################################################
################################################################################################

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
sample_minusE2.bed.srm.dir/ \
sample_plusE2.bed.srm.dir/ \
> tss.hg38.counts.NORM1e7.size.1kb.minusE2_plusE2.txt

################################################################################################
################################################################################################