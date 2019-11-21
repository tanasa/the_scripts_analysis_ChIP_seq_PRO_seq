#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################

library("ggplot2")
library("reshape2")
library("edgeR")
library("plyr")
library("dplyr")

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
###################################################################### initially, a few settings :

SIZE=2
LIMITS=c(0,150)
COLORS=c("red", "blue")

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
### to read the file with ALL the GENES on hg38 with the RAW COUNTS, as computed by HOMER 

ALL = "tss.hg38.counts.RAW.txt"

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
### here uploading the SETS of GENES :
### i.e. the longest isoform for ALL the RefSeq GENES on hg38, 
### the UP genes with PR, 
### the UP genes with no PR, 
### a set of RANDOMLY selected genes  

genes.all = read.delim("set.genes.here.ALL.longest.isoforms.TSS.HG38.bed", header = FALSE, sep = "\t", stringsAsFactors=F)
dim(genes.all)
head(genes.all)

genes.PR = read.delim("set.genes.UP.PAUSE.RELEASE.TSS.HG38.bed", header = FALSE, sep = "\t", stringsAsFactors=F)
dim(genes.PR)
head(genes.PR)

genes.OG = read.delim("set.genes.UP.other.genes.TSS.HG38.bed", header = FALSE, sep = "\t", stringsAsFactors=F)
dim(genes.OG)
head(genes.OG)

genes.RAND2 = read.delim("set.genes.here.random837.TSS.HG38.bed", header = FALSE, sep = "\t", stringsAsFactors=F)
dim(genes.RAND2)
head(genes.RAND2)

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
### SELECTING the RefSeq ID or the GENE NAMES 

genes.all.simple   = data.frame(RefSeq=genes.all$V5) 
genes.PR.simple    = data.frame(RefSeq=genes.PR$V4)
genes.OG.simple    = data.frame(RefSeq=genes.OG$V4)
genes.RAND2.simple  = data.frame(RefSeq=genes.RAND2$V4)

dim(genes.all.simple) 
dim(genes.PR.simple)  
dim(genes.OG.simple)  
dim(genes.RAND2.simple) 
 
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
### READING the file with the RAW COUNTS :
#############################################################################################################

all <- read.delim(ALL, header=T, sep = "\t", stringsAsFactors = F)

dim(all)
colnames(all)

#############################################################################################################
#############################################################################################################
############################### selecting the samples / the columns that we are interested in :
############################### typically in the HOMER file, the COUNTS are on the columns 20, 21, 22, 23, ..
#############################################################################################################

x.small <- data.frame(RefSeq = all[,1],
                      gene = all[,16],
                      minus = all[,21],
                      plus = all[,23])    

dim(x.small)
head(x.small)

#############################################################################################################
#############################################################################################################
############################## we do have multiple isoforms per gene, we work with the longest isoform / gene 
#############################################################################################################

x.small.LI <- merge(genes.all.simple,
                    x.small,
                    by.x = "RefSeq",
                    by.y = "RefSeq",
                    # all.x = FALSE,
                    # all.y = FALSE,
                    all = FALSE )

dim(x.small.LI) 
head(x.small.LI)

#############################################################################################################
############################################################################################################# 
######################### setting the DATAFRAME for NORMALIZATION and working with edgeR :
#############################################################################################################

TITLE_HERE="about.ALL.genes.hg38"

png(paste(TITLE_HERE, "display.RAW.values.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(x.small.LI$minus,
        x.small.LI$plus,
        ylim=LIMITS,
        col="orange", 
        main="before normalization")
dev.off()

### preparing the dataframe for NORMALIZATION :

row.names(x.small.LI) <- x.small.LI$gene

x.small.LI <- x.small.LI[,-1]
head(x.small.LI)
x.small.LI <- x.small.LI[,-1]
head(x.small.LI)

#############################################################################################################
########################## working with edgeR
#############################################################################################################

group <- factor(c("minus", "plus"))

libSizes <- as.vector(colSums(x.small.LI))

y <- DGEList(counts=x.small.LI, group=group, lib.size=libSizes)

#############################################################################################################
#############################################################################################################

# y <- calcNormFactors(y)
# counts.per.m <- cpm(y, normalized.lib.sizes=TRUE)

#############################################################################################################
#############################################################################################################
### filtering the entries with very low counts : 

z <- y

keep <- rowSums(cpm(z) > 3) >= 1
z <- z[keep, ,keep.lib.sizes=FALSE]
dim(z)

z <- calcNormFactors(z)
counts.per.m.z <- cpm(z, normalized.lib.sizes=TRUE)
counts.per.m.z <- as.data.frame(counts.per.m.z)

#############################################################################################################
#############################################################################################################
###################### display the NORMALIZED DATA :

png(paste(TITLE_HERE, "display.NORMALIZED.values.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(counts.per.m.z[,1], 
        counts.per.m.z[,2],
        ylim=LIMITS,
        col="orange", 
        main="after normalization")
dev.off()

#############################################################################################################
#############################################################################################################
### writing the file with the NORMALIZED COUNTS :

write.table(counts.per.m.z, 
            file=paste(TITLE_HERE,"filtered.and.normalized.CPM.txt", sep="."),
            quote=F, sep="\t", row.names=TRUE)

#############################################################################################################
#############################################################################################################
### WILCOXON test :
wilcox.test(counts.per.m.z$minus, counts.per.m.z$plus)

##############################################################################################################
##############################################################################################################
##############################################################################################################
### processing these DATAFRAMES :

counts.per.m.z <- as.data.frame(counts.per.m.z)
counts.per.m.z$gene <- row.names(counts.per.m.z) 

head(counts.per.m.z)

##############################################################################################################
##############################################################################################################
##############################################################################################################
######################################################################T#######################################
##############################################################################################################
############################################################################################################## 
### using counts.per.m.z
### and to integrate with other dataframes :
### genes.PR.simple  
### genes.OG.simple  
### genes.RAND2.simple 

length(genes.PR.simple$RefSeq)  
length(genes.OG.simple$RefSeq)  
length(genes.RAND2.simple$RefSeq) 

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## genes.PR

TITLE_HERE="genes.UP.PR"
COLORS=c("red","blue")
LIMITS=c(0,150)

genes.PR.simple.NORM <- merge(genes.PR.simple,
                              counts.per.m.z,
                              by.x = "RefSeq",
                              by.y = "gene",
                              # all.x = FALSE,
                              # all.y = FALSE,
                              all = FALSE )

dim(genes.PR.simple) 
dim(counts.per.m.z) 
dim(genes.PR.simple.NORM)  
head(genes.PR.simple.NORM)

########################################################################################################
########################################################################################################

png(paste(TITLE_HERE, "after.normalization.genes.PR.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.NORM$minus,
        genes.PR.simple.NORM$plus,
        ylim=LIMITS,
        col=COLORS, 
        notch=TRUE,
        ylab="normalized counts",
        main="genes.UP.PR")
dev.off()

### writing the file with the normalized counts :
write.table(genes.PR.simple.NORM,
            file=paste(TITLE_HERE,"after.normalization.genes.PR.txt",sep="."), sep="\t", row.names=TRUE)

### wilcoxon.test :
wilcox.test(genes.PR.simple.NORM$minus,genes.PR.simple.NORM$plus)

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## genes.OG

TITLE_HERE="genes.UP.OG"
COLORS=c("red","blue")
LIMITS=c(0,150)

genes.OG.simple.NORM <- merge(genes.OG.simple,
                              counts.per.m.z,
                              by.x = "RefSeq",
                              by.y = "gene",
                              # all.x = FALSE,
                              # all.y = FALSE,
                              all = FALSE )

dim(genes.OG.simple) 
dim(counts.per.m.z) 
dim(genes.OG.simple.NORM)   
head(genes.OG.simple.NORM)

########################################################################################################
########################################################################################################

png(paste(TITLE_HERE, "after.normalization.genes.OG.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.OG.simple.NORM$minus,
        genes.OG.simple.NORM$plus,
        ylim=LIMITS,
        col=COLORS, 
        notch=TRUE,
        ylab="normalized counts",
        main="genes.UP.OG")
dev.off()

### writing the file with the normalized counts :

write.table(genes.OG.simple.NORM,
            file=paste(TITLE_HERE,"after.normalization.genes.OG.txt",sep="."), sep="\t", row.names=TRUE)

### wilcoxon.test :
wilcox.test(genes.OG.simple.NORM$minus, genes.OG.simple.NORM$plus)

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## genes.RANDOM

TITLE_HERE="genes.RANDOMLY.selected"
COLORS=c("red","blue")
LIMITS=c(0,150)

genes.RAND2.simple.NORM <- merge(genes.RAND2.simple,
                                 counts.per.m.z,
                                 by.x = "RefSeq",
                                 by.y = "gene",
                                 # all.x = FALSE,
                                 # all.y = FALSE,
                                 all = FALSE )

dim(genes.RAND2.simple) 
dim(counts.per.m.z)  
dim(genes.RAND2.simple.NORM)  
head(genes.RAND2.simple.NORM)

########################################################################################################
########################################################################################################

png(paste(TITLE_HERE, "after.normalization.genes.RANDOM.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.RAND2.simple.NORM$minus,
        genes.RAND2.simple.NORM$plus,
        ylim=LIMITS,
        col=COLORS, 
        notch=TRUE,
        ylab="normalized counts",
        main="genes.RANDOMLY.selected")
dev.off()

### writing the file with the normalized counts :
write.table(genes.RAND2.simple.NORM,
            file=paste(TITLE_HERE,"after.normalization.genes.RANDOM.txt",sep="."), sep="\t", row.names=TRUE)

### wilcoxon.test :
wilcox.test(genes.RAND2.simple.NORM$minus, genes.RAND2.simple.NORM$plus)

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## BOXPLOTS of 4

TITLE_HERE="comp:genes.PR:genes.non.PR"
COLORS2 = c("red", "blue", "red", "blue")

pdf(paste(TITLE_HERE, "display.the.BOXPLOTS.of.CPM", "pdf", sep="."), height=4.72, width=3.14)
boxplot(genes.PR.simple.NORM$minus, 
        genes.PR.simple.NORM$plus,
        genes.OG.simple.NORM$minus, 
        genes.OG.simple.NORM$plus,
        ylim=LIMITS, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="normalized counts")
dev.off()  

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.CPM", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.NORM$minus, 
        genes.PR.simple.NORM$plus,
        genes.OG.simple.NORM$minus, 
        genes.OG.simple.NORM$plus,
        ylim=LIMITS, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="normalized counts")
dev.off() 


median(genes.PR.simple.NORM$minus) 
median(genes.PR.simple.NORM$plus)
median(genes.OG.simple.NORM$minus) 
median(genes.OG.simple.NORM$plus)

wilcox.test(genes.PR.simple.NORM$minus, genes.PR.simple.NORM$plus) 
wilcox.test(genes.OG.simple.NORM$minus, genes.OG.simple.NORM$plus) 
wilcox.test(genes.PR.simple.NORM$minus, genes.OG.simple.NORM$minus) 
wilcox.test(genes.PR.simple.NORM$plus, genes.OG.simple.NORM$plus) 

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

TITLE_HERE="comp:genes.PR:genes.RANDOM"
COLORS2 = c("red", "blue", "red", "blue")

pdf(paste(TITLE_HERE, "display.the.BOXPLOTS.of.CPM", "pdf", sep="."), height=4.72, width=3.14)
boxplot(genes.PR.simple.NORM$minus, 
        genes.PR.simple.NORM$plus,
        genes.RAND2.simple.NORM$minus, 
        genes.RAND2.simple.NORM$plus,
        ylim=LIMITS, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="normalized counts")
dev.off()  

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.CPM", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.NORM$minus, 
        genes.PR.simple.NORM$plus,
        genes.RAND2.simple.NORM$minus, 
        genes.RAND2.simple.NORM$plus,
        ylim=LIMITS, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="normalized counts")
dev.off()  

 
median(genes.PR.simple.NORM$minus) 
median(genes.PR.simple.NORM$plus)
median(genes.RAND2.simple.NORM$minus) 
median(genes.RAND2.simple.NORM$plus)

wilcox.test(genes.PR.simple.NORM$minus, genes.PR.simple.NORM$plus)
wilcox.test(genes.RAND2.simple.NORM$minus, genes.RAND2.simple.NORM$plus)
wilcox.test(genes.PR.simple.NORM$minus, genes.RAND2.simple.NORM$minus)
wilcox.test(genes.PR.simple.NORM$plus, genes.RAND2.simple.NORM$plus)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
##### here printing the MEDIAN VALUES and the P_VALUES in another FILE :

TITLE_HERE="comp:here:the.MEDIANS:and:the.PVALUES.txt"

#########################################################################################################
############ the comparison UP_PR vs UP_nonPR
#########################################################################################################

write.table("\n UP_PR :vs: UP_nonPR (OG) \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("median_ChIP:genes_UP_PR:minus", median(genes.PR.simple.NORM$minus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("median_ChIP:genes_UP_PR:plus", median(genes.PR.simple.NORM$plus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("median_ChIP:genes_OG:minus", median(genes.OG.simple.NORM$minus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("median_ChIP:genes_OG:plus", median(genes.OG.simple.NORM$plus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

#########################################################################################################

write.table("\n UP_PR :vs: UP_nonPR (OG) \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("wilcoxon_pval:genes_PR:minus:vs:genes_PR:plus", 
            wilcox.test(genes.PR.simple.NORM$minus, genes.PR.simple.NORM$plus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("wilcoxon_pval:genes_OG:minus:vs:genes_OG:plus", 
            wilcox.test(genes.OG.simple.NORM$minus, genes.OG.simple.NORM$plus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("wilcoxon_pval:genes_PR:minus:vs:genes_OG:minus", 
            wilcox.test(genes.PR.simple.NORM$minus, genes.OG.simple.NORM$minus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("wilcoxon_pval:genes_PR:plus:vs:genes_OG:plus", 
            wilcox.test(genes.PR.simple.NORM$plus, genes.OG.simple.NORM$plus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

#########################################################################################################
############ the comparison UP_PR vs RANDOM
#########################################################################################################

write.table("\n UP_PR :vs: RANDOM \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("median_ChIP:genes_UP_PR:minus", median(genes.PR.simple.NORM$minus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("median_ChIP:genes_UP_PR:plus", median(genes.PR.simple.NORM$plus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("median_ChIP:genes_RAND:minus", median(genes.RAND2.simple.NORM$minus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("median_ChIP:genes_RAND:plus", median(genes.RAND2.simple.NORM$plus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

#########################################################################################################

write.table("\n UP_PR :vs: RANDOM \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("wilcoxon_pval:genes_PR:minus:vs:genes_PR:plus", 
            wilcox.test(genes.PR.simple.NORM$minus, genes.PR.simple.NORM$plus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("wilcoxon_pval:genes_RAND:minus:vs:genes_RAND:plus", 
            wilcox.test(genes.RAND2.simple.NORM$minus, genes.RAND2.simple.NORM$plus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("wilcoxon_pval:genes_PR:minus:vs:genes_RAND:minus", 
            wilcox.test(genes.PR.simple.NORM$minus, genes.RAND2.simple.NORM$minus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)


write.table(paste("wilcoxon_pval:genes_PR:plus:vs:genes_RAND:plus", 
            wilcox.test(genes.PR.simple.NORM$plus, genes.RAND2.simple.NORM$plus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
