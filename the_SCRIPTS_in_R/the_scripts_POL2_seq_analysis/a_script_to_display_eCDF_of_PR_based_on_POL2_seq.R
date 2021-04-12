############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

library("ggplot2")
library("reshape2")
library("dplyr")
library("plyr")

######################################################################
######################################################################
######################################################################
######################################################################

breaks = c(0.05, 0.5, 5, 50)
LIMITS = c(0.05, 500)

######################################################################
COLUMNS_BEGIN = c(27,28)          #################################### it can be changed
###################################################################### there is a difference between the columns : 32-27 = 5
COLUMNS_LATE = c(32,33)           #################################### it can be changed
######################################################################

name_FILE_PVALUES="display.the.P_VALUES.KS.test.txt"    
SIZE = 2

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

ALL = "" ##################################### reading the INPUT FILE

all <- read.delim(ALL, header=TRUE, sep = "\t", stringsAsFactors = F)
colnames(all)

############################################################################################################################################
############################################################################################################################################ reading the GENES of hg38 
############################################################################################################################################

genes.all = read.delim("genes.all.TSS.HG38.bed",
                       header = FALSE, sep = "\t", stringsAsFactors=F)
dim(genes.all)
head(genes.all)

genes.PR = read.delim("genes.TSS.UP.PR.bed",
                      header = FALSE, sep = "\t", stringsAsFactors=F)
dim(genes.PR)
head(genes.PR)

genes.OG = read.delim("genes.TSS.UPnoPR.bed",
                      header = FALSE, sep = "\t", stringsAsFactors=F)
dim(genes.OG)
head(genes.OG)

genes.RAND2 = read.delim("genes.random837.bed",
                       header = FALSE, sep = "\t", stringsAsFactors=F)
dim(genes.RAND2)
head(genes.RAND2)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

TITLE_HERE="ALL.genes.hg38"

all$RefSeq <- all[,1]

all.small <- merge(genes.all,
                   all,
                   by.x = "V5",
                   by.y = "RefSeq",
                   all = FALSE )

dim(all.small) 
head(all.small)
colnames(all.small)
 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ 
############################################################################################################################################
############################################################################################################################################ 

all_PR <- all.small[, COLUMNS_BEGIN]      ################################################################################################## here COLUMNS_BEGIN      
colnames(all_PR) <- c("minus", "plus")

all_m <- melt(all_PR)
head(all_m)

ggplot(all_m, aes(value, colour = variable)) + stat_ecdf(size=SIZE) +   
                                               # scale_x_log10(breaks=breaks) +
                                               scale_x_log10(breaks=breaks, limits=LIMITS) +  
                                               theme(axis.text.x = element_text(size=8, angle=45)) + 
                                               labs(title=TITLE_HERE, y ="% genes", x="PR") +
                                               theme_classic()  +
                                               theme_bw() + 
            scale_colour_manual(values = c("red", "blue"), name="") +
            theme (axis.text = element_text(size=12),
            axis.title.x = element_text(size=12, vjust= -0.2),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.y = element_text(size=12, angle=90, vjust= NULL),
            legend.position = c(.80,0.18),
            legend.background = element_rect(),
            legend.title = element_blank(), 
            legend.key = element_blank(),
            legend.text = element_text(size = 12))

ggsave(paste("display.eCDF", TITLE_HERE, "genes", "png", sep=".")) 

ks.test(all_PR[,1], all_PR[,2])    

write.table(paste("KS p-value : ", "ALL genes", sep=""),  
            file=name_FILE_PVALUES, 
            sep="\t", append=T, quote=F, col.names=F, row.names=F)

write.table(ks.test(all_PR[,1], all_PR[,2])$p.value, 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F,col.names=F, row.names=F)

write.table("\n", 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F, col.names=F, row.names=F)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ genes.PR.simple.fused 

TITLE_HERE="genes.UP.with.PR"  

genes.PR.simple.fused <- merge(genes.PR,
                               all.small,
                               by.x = "V5",
                               by.y = "V5",
                               # all.x = FALSE,
                               # all.y = FALSE,
                               all = FALSE )

dim(genes.PR.simple.fused)
head(genes.PR.simple.fused)

genes.PR.simple.fused.PR <- genes.PR.simple.fused[, COLUMNS_LATE]  ######################################################################## here to change the COLUMNS
colnames(genes.PR.simple.fused.PR) <- c("minus", "plus")

write.table(genes.PR.simple.fused.PR, file=paste(ALL, "genes.PR.simple.fused", "txt", sep="."), 
                                      sep="\t", quote=F, row.names=F, col.names=F)        

all_m <- melt(genes.PR.simple.fused.PR)
head(all_m)

ggplot(all_m, aes(value, colour = variable)) + stat_ecdf(size=SIZE) +   
                                               # scale_x_log10(breaks=breaks) +
                                               scale_x_log10(breaks=breaks, limits=LIMITS) +  
                                               theme(axis.text.x = element_text(size=8, angle=45)) + 
                                               labs(title=TITLE_HERE, y ="% genes", x="PR") +
                                               theme_classic() +
                                               theme_bw() + 
            scale_colour_manual(values = c("red", "blue"), name="") +
            theme (axis.text = element_text(size=12),
            axis.title.x = element_text(size=12, vjust= -0.2),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.y = element_text(size=12, angle=90, vjust= NULL),
            legend.position = c(.80,0.18),
            legend.background = element_rect(),
            legend.title = element_blank(), 
            legend.key = element_blank(),
            legend.text = element_text(size = 12))

ggsave(paste("display.eCDF", TITLE_HERE, "genes", "png", sep="."))

ks.test(genes.PR.simple.fused.PR[,1], genes.PR.simple.fused.PR[,2]) 

write.table(paste("KS p-value : ", TITLE_HERE, sep=""),  
            file=name_FILE_PVALUES, 
            sep="\t", append=T, quote=F, col.names=F, row.names=F)
write.table(ks.test(genes.PR.simple.fused.PR[,1], genes.PR.simple.fused.PR[,2])$p.value, 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F,col.names=F, row.names=F)
write.table("\n", 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F, col.names=F, row.names=F)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ genes.OG.simple.fused

TITLE_HERE="genes.UP.Other.Genes" 

genes.OG.simple.fused <- merge(genes.OG,
                               all.small,
                               by.x = "V5",
                               by.y = "V5",
                               all = FALSE )

dim(genes.OG.simple.fused)
head(genes.OG.simple.fused)

genes.OG.simple.fused.PR <- genes.OG.simple.fused[, COLUMNS_LATE]   ####################################################################### here to change the COLUMNS 
colnames(genes.OG.simple.fused.PR) <- c("minus", "plus")

write.table(genes.OG.simple.fused.PR, 
            file=paste(ALL, "genes.OG.simple.fused", "txt", sep="."), 
            sep="\t", quote=F, row.names=F, col.names=F)        

all_m <- melt(genes.OG.simple.fused.PR)
head(all_m)

ggplot(all_m, aes(value, colour = variable)) + stat_ecdf(size=SIZE) +   
                                               # scale_x_log10(breaks=breaks) +
                                               scale_x_log10(breaks=breaks, limits=LIMITS) +  
                                               theme(axis.text.x = element_text(size=8, angle=45)) + 
                                               labs(title=TITLE_HERE, y ="% genes", x="PR") +
                                               theme_classic() +
                                               theme_bw() + 
            scale_colour_manual(values = c("red", "blue"), name="") +
            theme (axis.text = element_text(size=12),
            axis.title.x = element_text(size=12, vjust= -0.2),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.y = element_text(size=12, angle=90, vjust= NULL),
            legend.position = c(.80,0.18),
            legend.background = element_rect(),
            legend.title = element_blank(), 
            legend.key = element_blank(),
            legend.text = element_text(size = 12))

ggsave(paste("display.eCDF", TITLE_HERE, "genes", "png", sep="."))
 
ks.test(genes.OG.simple.fused.PR[,1], genes.OG.simple.fused.PR[,2])   # 

write.table(paste("KS p-value : ", TITLE_HERE, sep=""),  
            file=name_FILE_PVALUES, 
            sep="\t", append=T, quote=F, col.names=F, row.names=F)
write.table(ks.test(genes.OG.simple.fused.PR[,1], genes.OG.simple.fused.PR[,2])$p.value, 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F,col.names=F, row.names=F)
write.table("\n", 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F, col.names=F, row.names=F)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ genes.RAND2.simple.fused

TITLE_HERE="genes.RANDOM.set2" 

genes.RAND2.simple.fused <- merge(genes.RAND2,
                               all.small,
                               by.x = "V5",
                               by.y = "V5",
                               all = FALSE )

dim(genes.RAND2.simple.fused)
head(genes.RAND2.simple.fused)

genes.RAND2.simple.fused.PR <- genes.RAND2.simple.fused[, COLUMNS_LATE] ################################################################### here to change the COLUMNS
colnames(genes.RAND2.simple.fused.PR) <- c("minus", "plus")

write.table(genes.RAND2.simple.fused.PR, 
            file=paste(ALL, "genes.RAND2.simple.fused", "txt", sep="."), 
            sep="\t", quote=F, row.names=F, col.names=F)        

all_m <- melt(genes.RAND2.simple.fused.PR)
head(all_m)

ggplot(all_m, aes(value, colour = variable)) + stat_ecdf(size=SIZE) +   
                                               # scale_x_log10(breaks=breaks) +
                                               scale_x_log10(breaks=breaks, limits=LIMITS) +  
                                               theme(axis.text.x = element_text(size=8, angle=45)) + 
                                               labs(title=TITLE_HERE, y ="% genes", x="PR") +
                                               theme_classic() +
                                               theme_bw() + 
            scale_colour_manual(values = c("red", "blue"), name="") +
            theme (axis.text = element_text(size=12),
            axis.title.x = element_text(size=12, vjust= -0.2),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.y = element_text(size=12, angle=90, vjust= NULL),
            legend.position = c(.80,0.18),
            legend.background = element_rect(),
            legend.title = element_blank(), 
            legend.key = element_blank(),
            legend.text = element_text(size = 12))

ggsave(paste("display.eCDF", TITLE_HERE, "genes", "png", sep="."))

ks.test(genes.RAND2.simple.fused.PR[,1], genes.RAND2.simple.fused.PR[,2])   # 

write.table(paste("KS p-value : ", TITLE_HERE, sep=""),  
            file=name_FILE_PVALUES, 
            sep="\t", append=T, quote=F, col.names=F, row.names=F)
write.table(ks.test(genes.RAND2.simple.fused.PR[,1], genes.RAND2.simple.fused.PR[,2])$p.value, 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F,col.names=F, row.names=F)
write.table("\n", 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F, col.names=F, row.names=F)

###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
