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

TRIMMOMATIC="/home/btanasa/SOFTWARE/Trimmomatic-0.38"

################################################################################################

hg38="/home/btanasa/SOFTWARE/bowtie2_genomes/hg38/hg38_genome.bowtie2"

################################################################################################

FILE=$1

OUTPUT="${FILE}.after.trimmomatic.fastq.gz"   #### the output after trimming the reads

################################################################################################
################################################################################################
##### using FASTQC BEFORE TRIMMING 

mkdir "${FILE}.report.fastqc.before.trimming"

fastqc -t 12 \
-o "${FILE}.report.fastqc.before.trimming" \
$FILE 

################################################################################################
################################################################################################

java -jar $TRIMMOMATIC/trimmomatic-0.38.jar SE \
-phred33 -threads 12 \
$FILE \
$OUTPUT \
ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq_adapters_NEXTERA_from_BBMAP.txt:2:30:10 \
SLIDINGWINDOW:4:15  \
LEADING:10 \
TRAILING:10 \
MINLEN:36

################################################################################################
################################################################################################
##### using FASTQC AFTER TRIMMING 

mkdir "${FILE}.report.fastqc.here.trimmed"

fastqc -t 12 \
-o "${FILE}.report.fastqc.here.trimmed" \
$OUTPUT

###############################################################################################
###############################################################################################
##### the file that is the output of the trimming process :

FILE=$OUTPUT

##############################################################################################
###############################################################################################
##### the alignment with BOWTIE2 :

bowtie2 --very-sensitive -q -p 12 -x $hg38 -U $FILE -S "${FILE}.SAM"

###############################################################################################
###############################################################################################
##### generating a BED file from SAM :

samtools view "${FILE}.SAM" -Sb | bamToBed -i stdin > "${FILE}.bed"

###############################################################################################
###############################################################################################
##### making a TAG DIRECTORY and a UCSC file in HOMER, using all the aligned reads  

cut -f1,2,3,6 "${FILE}.bed" > "${FILE}.bed4"

makeTagDirectory "${FILE}.bed.dir" "${FILE}.bed4"
makeUCSCfile     "${FILE}.bed.dir" -o auto

###############################################################################################
###############################################################################################
##### keeping only 1_copy_read/genome_coordinate
 
sort -k1,1 -k2,2n -k3,3n -u "${FILE}.bed" > "${FILE}.bed.srm"

###############################################################################################
###############################################################################################
##### making a TAG DIRECTORY and a UCSC file in HOMER, 
##### using the new file in HOMER that has 1_copy_read/genome_coordinate

cut -f1,2,3,6 "${FILE}.bed.srm" > "${FILE}.bed.srm4"

makeTagDirectory "${FILE}.bed.srm.dir" "${FILE}.bed.srm4" 
makeUCSCfile     "${FILE}.bed.srm.dir" -o auto

###############################################################################################
###############################################################################################
##### removing some files that may not be needed

rm "${FILE}.SAM"
rm "${FILE}.bed4"
rm "${FILE}.bed.srm4" 

###############################################################################################
###############################################################################################
