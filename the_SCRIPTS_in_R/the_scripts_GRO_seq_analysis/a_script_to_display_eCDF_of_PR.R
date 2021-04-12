############################################################################################################+
#############################################################################################################
#############################################################################################################
#############################################################################################################

library("ggplot2")
library("reshape2")
library("dplyr")
library("plyr")

############################################################################################################+
#############################################################################################################
#############################################################################################################
#############################################################################################################

breaks = c(0.05, 0.5, 5, 50)
LIMITS = c(0.05, 500)
SIZE = 2

######################################################################
######################################################################

ALL = "" ################################### the file with info on ALL the GENES of hg18 (RPKM, PAUSING, RAW) 

ALL.info.hg18 = read.delim(ALL, header = TRUE, sep = "\t", stringsAsFactors=F)

dim(ALL.info.hg18)
head(ALL.info.hg18)
colnames(ALL.info.hg18)

######################################################################
###################################################################### reading the GENES of hg38 annotations

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

############################################################################################################+
#############################################################################################################
#############################################################################################################
#############################################################################################################

TITLE_HERE="genes.ALL.genes.hg18"

all <- read.delim(ALL, header=T, sep = "\t", stringsAsFactors = F)

all_PR <- all[, c(51,54)]
colnames(all_PR) <- c("minus", "plus")

ks.test(all_PR[,1], all_PR[,2])

############################################################################################################+
#############################################################################################################
#############################################################################################################
#############################################################################################################

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

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

genes.all.simple    = data.frame(RefSeq=genes.all$V4)
genes.PR.simple     = data.frame(RefSeq=genes.PR$V4)
genes.OG.simple     = data.frame(RefSeq=genes.OG$V4)
genes.RAND2.simple  = data.frame(RefSeq=genes.RAND2$V4)

dim(genes.all.simple) 
dim(genes.PR.simple)  
dim(genes.OG.simple)  
dim(genes.RAND2.simple) 

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## genes.all.simple

TITLE_HERE="genes.ALL.hg38.hg18"  

genes.all.simple.fused <- merge(genes.all.simple,
                                ALL.info.hg18,
                                by.x = "RefSeq",
                                by.y = "Symbol.x",
                                # all.x = FALSE,
                                # all.y = FALSE,
                                all = FALSE )

dim(genes.all.simple.fused)
head(genes.all.simple.fused)

genes.all.simple.fused.PR <- genes.all.simple.fused[, c(51,54)]
colnames(genes.all.simple.fused.PR) <- c("minus", "plus")

all_m <- melt(genes.all.simple.fused.PR)
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

ks.test(genes.all.simple.fused.PR[,1], genes.all.simple.fused.PR[,2], paired = TRUE)$p.value
                                                                                                                                                                                                                    
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################                                                                                                                                                                                                        
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## genes.PR

TITLE_HERE="genes.UP.with.PR"  

genes.PR.simple.fused <- merge(genes.PR.simple,
                               ALL.info.hg18,
                               by.x = "RefSeq",
                               by.y = "Symbol.x",
                               # all.x = FALSE,
                               # all.y = FALSE,
                               all = FALSE )

dim(genes.PR.simple.fused)
head(genes.PR.simple.fused)

genes.PR.simple.fused.PR <- genes.PR.simple.fused[, c(51,54)]
colnames(genes.PR.simple.fused.PR) <- c("minus", "plus")

all_m <- melt(genes.PR.simple.fused.PR)
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

ks.test(genes.PR.simple.fused.PR[,1], genes.PR.simple.fused.PR[,2], paired = TRUE)$p.value

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## genes.OG

TITLE_HERE="genes.UP.other.nonPR" 

genes.OG.simple.fused <- merge(genes.OG.simple,
                               ALL.info.hg18,
                               by.x = "RefSeq",
                               by.y = "Symbol.x",
                               # all.x = FALSE,
                               # all.y = FALSE,
                               all = FALSE )

dim(genes.OG.simple.fused)
head(genes.OG.simple.fused)

genes.OG.simple.fused.PR <- genes.OG.simple.fused[, c(51,54)]
colnames(genes.OG.simple.fused.PR) <- c("minus", "plus")

all_m <- melt(genes.OG.simple.fused.PR)
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

ks.test(genes.OG.simple.fused.PR[,1], genes.OG.simple.fused.PR[,2], paired = TRUE)$p.value

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## genes.RAND2

TITLE_HERE="genes.RANDOM.set2" 

genes.RAND2.simple.fused <- merge(genes.RAND2.simple,
                               ALL.info.hg18,
                               by.x = "RefSeq",
                               by.y = "Symbol.x",
                               # all.x = FALSE,
                               # all.y = FALSE,
                               all = FALSE )

dim(genes.RAND2.simple.fused)
head(genes.RAND2.simple.fused)

genes.RAND2.simple.fused.PR <- genes.RAND2.simple.fused[, c(51,54)]
colnames(genes.RAND2.simple.fused.PR) <- c("minus", "plus")

all_m <- melt(genes.RAND2.simple.fused.PR)
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

ks.test(genes.RAND2.simple.fused.PR[,1], genes.RAND2.simple.fused.PR[,2], paired = TRUE)$p.value   

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
######################################################################################################### the KS tests

write.table("the RESULTS of KS TESTS : ALL genes", 
            file = "the.results.KS.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(ks.test(genes.all.simple.fused.PR[,1], genes.all.simple.fused.PR[,2], paired = TRUE)$p.value, 
            file = "the.results.KS.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table("the RESULTS of KS TESTS : UP_PR", 
            file = "the.results.KS.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(ks.test(genes.PR.simple.fused.PR[,1], genes.PR.simple.fused.PR[,2], paired = TRUE)$p.value, 
            file = "the.results.KS.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table("the RESULTS of KS TESTS : UPnonPR", 
            file = "the.results.KS.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(ks.test(genes.OG.simple.fused.PR[,1], genes.OG.simple.fused.PR[,2], paired = TRUE)$p.value, 
            file = "the.results.KS.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table("the RESULTS of KS TESTS : RAND2", 
            file = "the.results.KS.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(ks.test(genes.RAND2.simple.fused.PR[,1], genes.RAND2.simple.fused.PR[,2], paired = TRUE)$p.value, 
            file = "the.results.KS.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
