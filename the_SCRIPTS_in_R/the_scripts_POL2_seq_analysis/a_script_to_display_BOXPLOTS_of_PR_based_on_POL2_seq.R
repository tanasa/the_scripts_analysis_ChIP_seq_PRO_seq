############################################################################################################################################
############################################################################################################################################

library("ggplot2")
library("reshape2")
library("dplyr")
library("plyr")

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

breaks = c(0.05, 0.5, 5, 50)
LIMITS = c(0.05, 500)

############################################################ to note a difference 20 - 15 = 5 between the COLUMNS in the file

COLUMNS_BEGIN = c(27,28)                                #### to change the columns
COLUMNS_LATE = c(32,33)                                 #### to change the columns 

name_FILE_PVALUES="display.the.P_VALUES.WILCOXON.test.txt"

############################################################################################################################################
############################################################################################################################################

ALL = "" ###################### reading the file of POL2 data on hg38

all <- read.delim(ALL, header=TRUE, sep = "\t", stringsAsFactors = F)
colnames(all)

############################################################################################################################################
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

SIZE = 2
LIMITS = c(0,15)
COLORS = c("red", "blue")

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

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
############################################################################################################################################

TITLE_HERE="ALL.genes.hg38"

all_PR <- all.small[, COLUMNS_BEGIN]            ############################################################################################ here to change THE COLUMNS
colnames(all_PR) <- c("minus", "plus")

############################################################################### here the BOXPLOTS in R

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(all_PR$minus, 
        all_PR$plus,
        ylim=LIMITS, 
        col=COLORS, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off()          

############################################################################################
wilcox.test(all_PR[,1], all_PR[,2])$p.value
############################################################################################

write.table(paste("WILCOXON p-value : ", TITLE_HERE, sep=""),  
            file=name_FILE_PVALUES, 
            sep="\t", append=T, quote=F, col.names=F, row.names=F)
write.table(wilcox.test(all_PR[,1], all_PR[,2])$p.value, 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F,col.names=F, row.names=F)
write.table("\n", 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F, col.names=F, row.names=F)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ genes.PR.simple.fused
                                                                                                                                            
TITLE_HERE="genes.UP.with.PR"  

genes.PR.simple.fused <- merge(genes.PR,
                               all.small,
                               by.x = "V5",
                               by.y = "V5",
                               all = FALSE )

dim(genes.PR.simple.fused)
# head(genes.PR.simple.fused)

genes.PR.simple.fused.PR <- genes.PR.simple.fused[, COLUMNS_LATE]  ####################################################################### here to change the COLUMNS
colnames(genes.PR.simple.fused.PR) <- c("minus", "plus")

############################################################################### here the BOXPLOTS in R

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=6, units="cm", res=800)
boxplot(genes.PR.simple.fused.PR$minus, 
        genes.PR.simple.fused.PR$plus,
        ylim=LIMITS, 
        col=COLORS, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off()          

############################################################################################
wilcox.test(genes.PR.simple.fused.PR[,1], genes.PR.simple.fused.PR[,2])$p.value
############################################################################################

write.table(paste("WILCOXON p-value : ", TITLE_HERE, sep=""),  
            file=name_FILE_PVALUES, 
            sep="\t", append=T, quote=F, col.names=F, row.names=F)
write.table(wilcox.test(genes.PR.simple.fused.PR[,1], genes.PR.simple.fused.PR[,2])$p.value, 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F,col.names=F, row.names=F)
write.table("\n", 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F, col.names=F, row.names=F)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ genes.OG.simple.fused 

TITLE_HERE="genes.UP.with.Other.Genes" 

genes.OG.simple.fused <- merge(genes.OG,
                               all.small,
                               by.x = "V5",
                               by.y = "V5",
                               all = FALSE )

dim(genes.OG.simple.fused)
# head(genes.OG.simple.fused)

genes.OG.simple.fused.PR <- genes.OG.simple.fused[, COLUMNS_LATE]  ####################################################################### here to change the COLUMNS 
colnames(genes.OG.simple.fused.PR) <- c("minus", "plus")

############################################################################### here the BOXPLOTS in R

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.OG.simple.fused.PR$minus, 
        genes.OG.simple.fused.PR$plus,
        ylim=LIMITS, 
        col=COLORS, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off()   

############################################################################################
wilcox.test(genes.OG.simple.fused.PR[,1], genes.OG.simple.fused.PR[,2])$p.value
############################################################################################

write.table(paste("WILCOXON p-value : ", TITLE_HERE, sep=""),  
            file=name_FILE_PVALUES, 
            sep="\t", append=T, quote=F, col.names=F, row.names=F)
write.table(wilcox.test(genes.OG.simple.fused.PR[,1], genes.OG.simple.fused.PR[,2])$p.value, 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F,col.names=F, row.names=F)
write.table("\n", 
            file=name_FILE_PVALUES, 
            sep="\n", append=T, quote=F, col.names=F, row.names=F)

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

genes.RAND2.simple.fused.PR <- genes.RAND2.simple.fused[, COLUMNS_LATE] #################################################################### here to change the COLUMNS
colnames(genes.RAND2.simple.fused.PR) <- c("minus", "plus")

############################################################################### here the BOXPLOTS in R

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.RAND2.simple.fused.PR$minus, 
        genes.RAND2.simple.fused.PR$plus,
        ylim=LIMITS, 
        col=COLORS, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off()          

############################################################################################
wilcox.test(genes.RAND2.simple.fused.PR[,1], genes.RAND2.simple.fused.PR[,2])$p.value
############################################################################################

write.table(paste("WILCOXON p-value : ", TITLE_HERE, sep=""),  
            file=name_FILE_PVALUES, 
            sep="\t", append=T, quote=F, col.names=F, row.names=F)
write.table(wilcox.test(genes.RAND2.simple.fused.PR[,1], genes.RAND2.simple.fused.PR[,2])$p.value, 
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
############################################################################################################################################ to display the BOXPLOTS of 4  

TITLE_HERE="genes.PR:genes.non.PR"
COLORS2 = c("red", "blue", "red", "blue")
LIMITS2 = c(0,15)

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.fused.PR$minus, 
        genes.PR.simple.fused.PR$plus,
        genes.OG.simple.fused.PR$minus, 
        genes.OG.simple.fused.PR$plus,
        ylim=LIMITS2, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off() 
 
##################################################################

# write.table(paste("WILCOXON p-value : ", TITLE_HERE, sep=""),  
#            file=name_FILE_PVALUES, 
#             sep="\t", append=T, quote=F, col.names=F, row.names=F)
# write.table("1st and 2nd :", 
#             file=name_FILE_PVALUES, 
#             sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table(wilcox.test(genes.PR.simple.fused.PR$minus, genes.PR.simple.fused.PR$plus)$p.value, 
#             file=name_FILE_PVALUES, 
#             sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table("3rd and 4th :", 
#             file=name_FILE_PVALUES, 
#             sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table(wilcox.test(genes.OG.simple.fused.PR$minus, genes.OG.simple.fused.PR$plus)$p.value, 
#             file=name_FILE_PVALUES, 
#             sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table("1st and 3rd :", 
#             file=name_FILE_PVALUES, 
#             sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table(wilcox.test(genes.PR.simple.fused.PR$minus, genes.OG.simple.fused.PR$minus)$p.value, 
#             file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table("2nd and 4th :", 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table(wilcox.test(genes.PR.simple.fused.PR$plus, genes.OG.simple.fused.PR$plus)$p.value, 
#             file=name_FILE_PVALUES, 
#             sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table("\n", 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ to display the BOXPLOTS of 4 

TITLE_HERE="genes.PR:genes.RAND"
COLORS2 = c("red", "blue", "red", "blue")
LIMITS2 = c(0,15)

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.fused.PR$minus, 
        genes.PR.simple.fused.PR$plus,
        genes.RAND2.simple.fused.PR$minus, 
        genes.RAND2.simple.fused.PR$plus,
        ylim=LIMITS2, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off()  
 
##################################################################

# write.table(paste("WILCOXON p-value : ", TITLE_HERE, sep=""),  
#            file=name_FILE_PVALUES, 
#            sep="\t", append=T, quote=F, col.names=F, row.names=F)
# write.table("1st and 2nd :", 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table(wilcox.test(genes.PR.simple.fused.PR$minus, genes.PR.simple.fused.PR$plus)$p.value, 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table("3rd and 4th :", 
#             file=name_FILE_PVALUES, 
#             sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table(wilcox.test(genes.RAND2.simple.fused.PR$minus, genes.RAND2.simple.fused.PR$plus)$p.value, 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table("1st and 3rd :", 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table(wilcox.test(genes.PR.simple.fused.PR$minus, genes.RAND2.simple.fused.PR$minus)$p.value, 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table("2nd and 4th :", 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table(wilcox.test(genes.PR.simple.fused.PR$plus, genes.RAND2.simple.fused.PR$plus)$p.value, 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)
# write.table("\n", 
#            file=name_FILE_PVALUES, 
#            sep="\n", append=T, quote=F,col.names=F, row.names=F)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
