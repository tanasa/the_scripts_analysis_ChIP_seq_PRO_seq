#############################################################################################################
#############################################################################################################
#############################################################################################################
############################################################################################################# 

library("ggplot2")
library("reshape2")
library("edgeR")
library("plyr")
library("dplyr")

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

ALL = ""  ### reading the INPUT FILE 

######################################################################
######################################################################

SIZE=2
LIMITS=c(0,400)
COLORS=c("red", "blue")

######################################################################
###################################################################### reading the GENES of hg38 

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

#############################################################################################################
#############################################################################################################
#############################################################################################################
############################################################################################################# 

genes.all.simple   = data.frame(RefSeq=genes.all$V5) ### we keep RefSeq
genes.PR.simple    = data.frame(RefSeq=genes.PR$V4)
genes.OG.simple    = data.frame(RefSeq=genes.OG$V4)
genes.RAND2.simple = data.frame(RefSeq=genes.RAND2$V4)

dim(genes.all.simple) 
dim(genes.PR.simple)  
dim(genes.OG.simple)  
dim(genes.RAND2.simple) 

############################################################################################################3
#############################################################################################################
#############################################################################################################
#############################################################################################################

TITLE_HERE = "ALL.genes.hg38"
all <- read.delim(ALL, header=T, sep = "\t", stringsAsFactors = F)

dim(all)
colnames(all)

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################ selecting the columns that we are interested in :
############################################################## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! very important - where the COLUMNS are changing !

x.small <- data.frame(RefSeq = all[,1],
                      gene = all[,9],
                      minus = all[,17],
                      plus = all[,18])    

dim(x.small)
head(x.small)

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

x.small.LI <- merge(genes.all.simple,
                    x.small,
                    by.x = "RefSeq",
                    by.y = "RefSeq",
                    all = FALSE )

dim(x.small.LI) 
head(x.small.LI)

############################################################################################################
############################################################################################################ 

row.names(x.small.LI) <- x.small.LI$gene

x.small.LI <- x.small.LI[,-1]
head(x.small.LI)
x.small.LI <- x.small.LI[,-1]
head(x.small.LI)

############################################################################################################
############################################################################################################ 
############################################################################################################ 

group <- factor(c("minus", "plus"))
libSizes <- as.vector(colSums(x.small.LI))

y <- DGEList(counts=x.small.LI, group=group, lib.size=libSizes)
z <- y

############################################################################################################
############################################################################################################

# y <- calcNormFactors(y)
# counts.per.m <- cpm(y, normalized.lib.sizes=TRUE)

############################################################################################################
############################################################################################################ FILTERING

keep <- rowSums(cpm(z) > 3) >= 1
z <- z[keep, ,keep.lib.sizes=FALSE]

dim(z)

############################################################################################################
############################################################################################################ NORMALIZATION 

z <- calcNormFactors(z)
counts.per.m.z <- cpm(z, normalized.lib.sizes=TRUE)
counts.per.m.z <- as.data.frame(counts.per.m.z)

###################### display the NORMALIZED DATA :

png(paste(TITLE_HERE, "ALL_genes", "display_after_NORMALIZATION.after.edgeR.filtering.png", sep="."))
boxplot(counts.per.m.z[,1], 
        counts.per.m.z[,2],
        ylim=LIMITS,
        col="yellow", 
        main="NORMALIZED COUNTS in edgeR after filtering")
dev.off()

### wilcoxon test :
wilcox.test(counts.per.m.z$minus, counts.per.m.z$plus)$p.value

##############################################################################################################
##############################################################################################################
##############################################################################################################
############################################################################################################## processing these DATAFRAMES :

counts.per.m.z <- as.data.frame(counts.per.m.z)
counts.per.m.z$gene <- row.names(counts.per.m.z) 

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
############################################################################################################## using counts.per.m.z

### the other dataframes are :

### genes.PR.simple  
### genes.OG.simple  
### genes.OC.simple 
### genes.RAND.simple 

length(counts.per.m.z$gene)
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
######################################################################################################## genes.PR.simple

TITLE_HERE="genes.PR"
COLORS=c("red","blue")
LIMITS=c(0,400)

genes.PR.simple.NORM <- merge(genes.PR.simple,
                              counts.per.m.z,
                              by.x = "RefSeq",
                              by.y = "gene",
                              all = FALSE )

dim(genes.PR.simple) 
dim(genes.PR.simple.NORM) 
dim(counts.per.m.z)    
head(genes.PR.simple.NORM)

########################################################################## printing the image into a PNG :

png(paste(TITLE_HERE, "counts.per.m.z.genes.PR.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.NORM$minus,
        genes.PR.simple.NORM$plus,
        ylim=LIMITS,
        col=COLORS, 
        notch=TRUE,
        ylab="normalized counts",
        main="genes.PR")
dev.off()

wilcox.test(genes.PR.simple.NORM$minus, genes.PR.simple.NORM$plus)$p.value

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## genes.OG.simple

TITLE_HERE="genes.OG"
COLORS=c("red","blue")
LIMITS=c(0,400)

genes.OG.simple.NORM <- merge(genes.OG.simple,
                              counts.per.m.z,
                              by.x = "RefSeq",
                              by.y = "gene",
                              # all.x = FALSE,
                              # all.y = FALSE,
                              all = FALSE )

dim(genes.OG.simple) 
dim(genes.OG.simple.NORM)   
dim(counts.per.m.z)  
head(genes.OG.simple.NORM)

########################################################################## printing the image into a PNG :

png(paste(TITLE_HERE, "counts.per.m.z.genes.OG.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.OG.simple.NORM$minus,
        genes.OG.simple.NORM$plus,
        ylim=LIMITS,
        col=COLORS, 
        notch=TRUE,
        ylab="normalized counts",
        main="genes.OG")
dev.off()

wilcox.test(genes.OG.simple.NORM$minus, genes.OG.simple.NORM$plus)$p.value

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
######################################################################################################## genes.RAND2.simple

TITLE_HERE="genes.RAND2"
COLORS=c("red","blue")
LIMITS=c(0,400)

genes.RAND2.simple.NORM <- merge(genes.RAND2.simple,
                              counts.per.m.z,
                              by.x = "RefSeq",
                              by.y = "gene",
                              all = FALSE )

dim(genes.RAND2.simple) 
dim(genes.RAND2.simple.NORM)  
dim(counts.per.m.z)    
head(genes.RAND2.simple.NORM)

########################################################################## printing the image into a PNG :

png(paste(TITLE_HERE, "counts.per.m.z.genes.RAND2.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.RAND2.simple.NORM$minus,
        genes.RAND2.simple.NORM$plus,
        ylim=LIMITS,
        col=COLORS, 
        notch=TRUE,
        ylab="normalized counts",
        main="genes.RAND2")
dev.off()

### wilcoxon.test :

wilcox.test(genes.RAND2.simple.NORM$minus, genes.RAND2.simple.NORM$plus)$p.value

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################### here to display the BOXPLOTS of 4  

TITLE_HERE="comp:genes.PR:genes.non.PR"
COLORS2 = c("red", "blue", "red", "blue")

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

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
############################################################################################################ here to display the BOXPLOTS of 4 

TITLE_HERE="comp:genes.PR:genes.RAND2"
COLORS2 = c("red", "blue", "red", "blue")

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

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
##### here printing the MEDIAN VALUES and the P_VALUES in a FILE :

TITLE_HERE="comp:here:the.MEDIANS:and:the.PVALUES.txt"

#############################################
############ the comparison UP_PR vs UP_nonPR

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

############

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

###################################################
############ the comparison UP_PR vs RANDOM

write.table("\n UP_PR :vs: RANDOM (r2) \n",
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

write.table(paste("median_ChIP:genes_RAND2:minus", median(genes.RAND2.simple.NORM$minus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("median_ChIP:genes_RAND2:plus", median(genes.RAND2.simple.NORM$plus), sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

#############

write.table("\n UP_PR :vs: RANDOM (r2) \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(paste("wilcoxon_pval:genes_PR:minus:vs:genes_PR:plus", 
            wilcox.test(genes.PR.simple.NORM$minus, genes.PR.simple.NORM$plus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("wilcoxon_pval:genes_RAND2:minus:vs:genes_RAND2:plus", 
            wilcox.test(genes.RAND2.simple.NORM$minus, genes.RAND2.simple.NORM$plus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("wilcoxon_pval:genes_PR:minus:vs:genes_RAND2:minus", 
            wilcox.test(genes.PR.simple.NORM$minus, genes.RAND2.simple.NORM$minus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

write.table(paste("wilcoxon_pval:genes_PR:plus:vs:genes_RAND2:plus", 
            wilcox.test(genes.PR.simple.NORM$plus, genes.RAND2.simple.NORM$plus)$p.value , sep=" : "),
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

####################################################################################################
####################################################################################################
