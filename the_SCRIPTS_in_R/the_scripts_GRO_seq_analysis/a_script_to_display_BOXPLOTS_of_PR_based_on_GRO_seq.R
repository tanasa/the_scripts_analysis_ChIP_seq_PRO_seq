############################################################################################################################################
############################################################################################################################################

library("ggplot2")
library("reshape2")
library("dplyr")
library("plyr")

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

SIZE = 2
LIMITS=c(0,20)
COLORS=c("red", "blue")

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

ALL = ""                                               #### a file with lots of info on ALL the GENES of hg18 (RAW, RPKM, PAUSING) 

ALL.info.hg18 = read.delim(ALL, header = TRUE, sep = "\t", stringsAsFactors=F)

dim(ALL.info.hg18)
head(ALL.info.hg18)
colnames(ALL.info.hg18)

############################################################################################################################################
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
############################################################################################################################################

TITLE_HERE="ALL.genes.hg18"
all <- read.delim(ALL, header=T, sep = "\t", stringsAsFactors = F)

all_PR <- all[, c(51,54)]                            ################################################################# in our case 51 and 54  
colnames(all_PR) <- c("minus", "plus")

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

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

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

genes.all.simple    = data.frame(RefSeq=genes.all$V4)
genes.PR.simple     = data.frame(RefSeq=genes.PR$V4)
genes.OG.simple     = data.frame(RefSeq=genes.OG$V4)
genes.RAND2.simple  = data.frame(RefSeq=genes.RAND2$V4)

dim(genes.all.simple) 
dim(genes.PR.simple)  
dim(genes.OG.simple)  
dim(genes.RAND2.simple) 
 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ genes.all.simple
############################################################################################################################################
############################################################################################################################################

TITLE_HERE="genes.ALL.hg38.hg18"  

genes.all.simple.fused <- merge(genes.all.simple,
                                ALL.info.hg18,
                                by.x = "RefSeq",
                                by.y = "Symbol.x",
                                all = FALSE )

dim(genes.all.simple.fused)
# head(genes.all.simple.fused)

genes.all.simple.fused.PR <- genes.all.simple.fused[, c(51,54)] 
colnames(genes.all.simple.fused.PR) <- c("minus", "plus")

wilcox.test(genes.all.simple.fused.PR$minus, genes.all.simple.fused.PR$plus)$p.value 

################################################################# displaying the BOXPLOTS (png) :

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.all.simple.fused.PR$minus, 
        genes.all.simple.fused.PR$plus,
        ylim=LIMITS, 
        col=COLORS, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off()          

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ genes.PR.simple.fused
############################################################################################################################################
############################################################################################################################################
                                                                                                                                                                                                                    
TITLE_HERE="genes.UP.with.PR"  

genes.PR.simple.fused <- merge(genes.PR.simple,
                               ALL.info.hg18,
                               by.x = "RefSeq",
                               by.y = "Symbol.x",
                               all = FALSE )

dim(genes.PR.simple.fused)
# head(genes.PR.simple.fused)

genes.PR.simple.fused.PR <- genes.PR.simple.fused[, c(51,54)]
colnames(genes.PR.simple.fused.PR) <- c("minus", "plus")

wilcox.test(genes.PR.simple.fused.PR$minus, genes.PR.simple.fused.PR$plus)$p.value  

################################################################# displaying the BOXPLOTS (png) :

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.fused.PR$minus, 
        genes.PR.simple.fused.PR$plus,
        ylim=LIMITS, 
        col=COLORS, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off()          

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ genes.OG.simple.fused
############################################################################################################################################
############################################################################################################################################

TITLE_HERE="genes.UP.Other.Genes" 

genes.OG.simple.fused <- merge(genes.OG.simple,
                               ALL.info.hg18,
                               by.x = "RefSeq",
                               by.y = "Symbol.x",
                               all = FALSE )

dim(genes.OG.simple.fused)
# head(genes.OG.simple.fused)

genes.OG.simple.fused.PR <- genes.OG.simple.fused[, c(51,54)]
colnames(genes.OG.simple.fused.PR) <- c("minus", "plus")

wilcox.test(genes.OG.simple.fused.PR$minus, genes.OG.simple.fused.PR$plus)$p.value 
  
################################################################# displaying the BOXPLOTS (png) :

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

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ genes.RAND2.simple.fused
############################################################################################################################################
############################################################################################################################################

TITLE_HERE="genes.RANDOM.set2" 

genes.RAND2.simple.fused <- merge(genes.RAND2.simple,
                               ALL.info.hg18,
                               by.x = "RefSeq",
                               by.y = "Symbol.x",
                               all = FALSE )

dim(genes.RAND2.simple.fused)
head(genes.RAND2.simple.fused)

genes.RAND2.simple.fused.PR <- genes.RAND2.simple.fused[, c(51,54)]
colnames(genes.RAND2.simple.fused.PR) <- c("minus", "plus")

wilcox.test(genes.RAND2.simple.fused.PR$minus, genes.RAND2.simple.fused.PR$plus)$p.value 
  
################################################################# displaying the BOXPLOTS (png) :

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

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ here to display the BOXPLOTS of 4
TITLE_HERE="genes.PR:genes.non.PR"
COLORS2 = c("red", "blue", "red", "blue")

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.fused.PR$minus, 
        genes.PR.simple.fused.PR$plus,
        genes.OG.simple.fused.PR$minus, 
        genes.OG.simple.fused.PR$plus,
        ylim=LIMITS, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off() 

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ here to display the BOXPLOTS of 4
TITLE_HERE="genes.PR:genes.RAND"
COLORS2 = c("red", "blue", "red", "blue")

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.PR", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.fused.PR$minus, 
        genes.PR.simple.fused.PR$plus,
        genes.RAND2.simple.fused.PR$minus, 
        genes.RAND2.simple.fused.PR$plus,
        ylim=LIMITS, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="PR")
dev.off()  
 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

write.table("the RESULTS of WILCOXON TESTS : ALL genes", 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(wilcox.test(genes.all.simple.fused.PR$minus, genes.all.simple.fused.PR$plus)$p.value, 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table("the RESULTS of WILCOXON TESTS : UP_PR", 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(wilcox.test(genes.PR.simple.fused.PR$minus, genes.PR.simple.fused.PR$plus)$p.value, 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table("the RESULTS of WILCOXON TESTS : UPnonPR", 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(wilcox.test(genes.OG.simple.fused.PR$minus, genes.OG.simple.fused.PR$plus)$p.value, 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table("the RESULTS of WILCOXON TESTS : RAND2", 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(wilcox.test(genes.RAND2.simple.fused.PR$minus, genes.RAND2.simple.fused.PR$plus)$p.value, 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
