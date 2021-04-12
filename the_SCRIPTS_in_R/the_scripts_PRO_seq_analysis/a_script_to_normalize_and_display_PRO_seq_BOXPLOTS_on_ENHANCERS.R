############################################################################################################################################
############################################################################################################################################

library("ggplot2")
library("reshape2")
library("edgeR")
library("plyr")
library("dplyr")

############################################################################################################################################
############################################################################################################################################
# A – 1248 ER bound enhancers, MEGA-TRANS, with induced eRNAs 
# B – 5694 ER bound, less active enhancers
# C – 11008 other enhancers

# 3.lists.enhancers.v3.bed
# 3.lists.enhancers.v3.bed.1248_ACTIVE.bed
# 3.lists.enhancers.v3.bed.LESS_ACTIVE_ENHANCERS.bed
# 3.lists.enhancers.v3.bed.OTHER_ENHANCERS.bed
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

ALL = ""   ############# reading the LIST of ENHANCERS with the RAW COUNTS

SIZE=2
LIMITS=c(0, 300)
COLORS=c("red", "blue")

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ reading the lists of ENHANCERS  
############################################################################################################################################
############################################################################################################################################

enhancers_1248_ACTIVE = read.delim("3.lists.enhancers.v3.bed.1248_ACTIVE.bed",
                                    header = FALSE, sep = "\t", stringsAsFactors=F)
dim(enhancers_1248_ACTIVE)
head(enhancers_1248_ACTIVE)

enhancers_LESS_ACTIVE = read.delim("3.lists.enhancers.v3.bed.LESS_ACTIVE_ENHANCERS.bed",
                                   header = FALSE, sep = "\t", stringsAsFactors=F)
dim(enhancers_LESS_ACTIVE)
head(enhancers_LESS_ACTIVE)

enhancers_OTHERS = read.delim("3.lists.enhancers.v3.bed.OTHER_ENHANCERS.bed",
                              header = FALSE, sep = "\t", stringsAsFactors=F)
dim(enhancers_OTHERS)
head(enhancers_OTHERS)

############################################################################################################################################
############################################################################################################################################

enhancers_1248_ACTIVE.simple    = data.frame(ENHANCERS = enhancers_1248_ACTIVE$V4)
enhancers_LESS_ACTIVE.simple    = data.frame(ENHANCERS = enhancers_LESS_ACTIVE$V4)
enhancers_OTHERS.simple    = data.frame(ENHANCERS = enhancers_OTHERS$V4)

dim(enhancers_1248_ACTIVE.simple) 
dim(enhancers_LESS_ACTIVE.simple)  
dim(enhancers_OTHERS.simple)  

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

TITLE_HERE="ALL.ENHANCERS.hg38"

all <- read.delim(ALL, header=T, sep = "\t", stringsAsFactors = F)
dim(all) 
colnames(all)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ selecting the COLUMNS

x.small <- data.frame(enhancer = all[,1],
                      gene = all[,16],
                      minus = all[,22],
                      plus = all[,23], 
                      minus1 = all[,24],
                      plus1 = all[,25]) 
 
dim(x.small)
head(x.small)

rownames(x.small) = x.small$enhancer
head(x.small)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

x.small.LI <- x.small 

dim(x.small.LI) 
head(x.small.LI)
 
x.small.LI <- x.small.LI[,-1]
head(x.small.LI)
x.small.LI <- x.small.LI[,-1]
head(x.small.LI)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ edgeR

group <- factor(c("minus", "plus", "minus1", "plus1"))
libSizes <- as.vector(colSums(x.small.LI))

y <- DGEList(counts=x.small.LI, group=group, lib.size=libSizes)
z <- y

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

# y <- calcNormFactors(y)
# counts.per.m <- cpm(y, normalized.lib.sizes=TRUE)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ FILTERING

keep <- rowSums(cpm(z) > 3) >= 1
z <- z[keep, ,keep.lib.sizes=FALSE]
dim(z)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ NORMALIZATION

z <- calcNormFactors(z)
counts.per.m.z <- cpm(z, normalized.lib.sizes=TRUE)
counts.per.m.z <- as.data.frame(counts.per.m.z)

###################### display the NORMALIZED DATA :

png(paste(TITLE_HERE, "ALL_ENHANCERS", "display_after_NORMALIZATION.after.edgeR.filtering.png", sep="."))
boxplot(counts.per.m.z[,1], 
        counts.per.m.z[,2],
        counts.per.m.z[,3],
        counts.per.m.z[,4],
        ylim=LIMITS,
        col="yellow", 
        main="NORMALIZED COUNTS in edgeR after filtering")
dev.off()

###################### wilcoxon test on the normalized counts :

wilcox.test(counts.per.m.z$minus1, counts.per.m.z$plus1)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

counts.per.m.z <- as.data.frame(counts.per.m.z)
counts.per.m.z$ENHANCERS <- row.names(counts.per.m.z) 
head(counts.per.m.z)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

dim(enhancers_1248_ACTIVE.simple) 
dim(enhancers_LESS_ACTIVE.simple)  
dim(enhancers_OTHERS.simple)  

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ ENHANCERS
############################################################################################################################################ 1248 ACTIVE

TITLE_HERE="enhancers.1248.ACTIVE"
COLORS=c("red","blue")
LIMITS=c(0,400)

enhancers_1248_ACTIVE.simple.NORM <- merge(enhancers_1248_ACTIVE.simple,
                                           counts.per.m.z,
                              by.x = "ENHANCERS",
                              by.y = "ENHANCERS",
                              all = FALSE )

dim(enhancers_1248_ACTIVE.simple) 
dim(enhancers_1248_ACTIVE.simple.NORM)  

###############################################################################################################################

png(paste(TITLE_HERE, "counts.per.m.z.enhancers_1248_ACTIVE.png", sep="."), height=12, width=7, units="cm", res=800)
boxplot(enhancers_1248_ACTIVE.simple.NORM$minus,
        enhancers_1248_ACTIVE.simple.NORM$plus,
        enhancers_1248_ACTIVE.simple.NORM$minus1,
        enhancers_1248_ACTIVE.simple.NORM$plus1,
        ylim=LIMITS,
        col=COLORS, 
        notch=TRUE,
        ylab="normalized counts",
        main="enhancers_1248_ACTIVE")
dev.off()

### wilcoxon.test :
wilcox.test(enhancers_1248_ACTIVE.simple.NORM$minus, enhancers_1248_ACTIVE.simple.NORM$plus)
wilcox.test(enhancers_1248_ACTIVE.simple.NORM$minus1, enhancers_1248_ACTIVE.simple.NORM$plus1)
wilcox.test(enhancers_1248_ACTIVE.simple.NORM$plus, enhancers_1248_ACTIVE.simple.NORM$plus1)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ ENHANCERS
############################################################################################################################################ LESS ACTIVE

TITLE_HERE="enhancers.LESS_ACTIVE"
COLORS=c("red","blue")
LIMITS=c(0,400)

enhancers_LESS_ACTIVE.simple.NORM <- merge(enhancers_LESS_ACTIVE.simple,
                                           counts.per.m.z,
                              by.x = "ENHANCERS",
                              by.y = "ENHANCERS",
                              all = FALSE )

dim(enhancers_LESS_ACTIVE.simple) 
dim(enhancers_LESS_ACTIVE.simple.NORM)  

################################################################################################################################

png(paste(TITLE_HERE, "counts.per.m.z.enhancers_LESS_ACTIVE.png", sep="."), height=12, width=7, units="cm", res=800)
boxplot(enhancers_LESS_ACTIVE.simple.NORM$minus,
        enhancers_LESS_ACTIVE.simple.NORM$plus,
        enhancers_LESS_ACTIVE.simple.NORM$minus1,
        enhancers_LESS_ACTIVE.simple.NORM$plus1,
        ylim=LIMITS,
        col=COLORS, 
        notch=TRUE,
        ylab="normalized counts",
        main="enhancers_LESS_ACTIVE")
dev.off()

### wilcoxon.test :
wilcox.test(enhancers_LESS_ACTIVE.simple.NORM$minus, enhancers_LESS_ACTIVE.simple.NORM$plus)
wilcox.test(enhancers_LESS_ACTIVE.simple.NORM$minus1, enhancers_LESS_ACTIVE.simple.NORM$plus1)
wilcox.test(enhancers_LESS_ACTIVE.simple.NORM$plus, enhancers_LESS_ACTIVE.simple.NORM$plus1)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ ENHANCERS
############################################################################################################################################ OTHERS

TITLE_HERE="enhancers.OTHERS"
COLORS=c("red","blue")
LIMITS=c(0,400)

enhancers_OTHERS.simple.NORM <- merge(enhancers_OTHERS.simple,
                                      counts.per.m.z,
                              by.x = "ENHANCERS",
                              by.y = "ENHANCERS",
                              all = FALSE )

dim(enhancers_OTHERS.simple) 
dim(enhancers_OTHERS.simple.NORM)  

################################################################################################################################

png(paste(TITLE_HERE, "counts.per.m.z.enhancers_OTHERS.png", sep="."), height=12, width=7, units="cm", res=800)
boxplot(enhancers_OTHERS.simple.NORM$minus,
        enhancers_OTHERS.simple.NORM$plus,
        enhancers_OTHERS.simple.NORM$minus1,
        enhancers_OTHERS.simple.NORM$plus1,
        ylim=LIMITS,
        col=COLORS, 
        notch=TRUE,
        ylab="normalized counts",
        main="enhancers_OTHERS")
dev.off()

### wilcoxon.test :
wilcox.test(enhancers_OTHERS.simple.NORM$minus, enhancers_OTHERS.simple.NORM$plus)
wilcox.test(enhancers_OTHERS.simple.NORM$minus1, enhancers_OTHERS.simple.NORM$plus1)
wilcox.test(enhancers_OTHERS.simple.NORM$plus, enhancers_OTHERS.simple.NORM$plus1)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ the BOXPLOTS of 4 

TITLE_HERE="comp:enhancers:1243:vs:LESS_ACTIVE.enhancers"
COLORS2 = c("red", "blue", "red", "blue")
LIMITS=c(0,400)

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.CPM", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(enhancers_1248_ACTIVE.simple.NORM$minus, 
        enhancers_1248_ACTIVE.simple.NORM$plus,
        enhancers_1248_ACTIVE.simple.NORM$minus1,
        enhancers_1248_ACTIVE.simple.NORM$plus1,

        enhancers_LESS_ACTIVE.simple.NORM$minus, 
        enhancers_LESS_ACTIVE.simple.NORM$plus,
        enhancers_LESS_ACTIVE.simple.NORM$minus1,
        enhancers_LESS_ACTIVE.simple.NORM$plus1,

        ylim=LIMITS, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="normalized counts")
dev.off() 

############################################################################################

median(enhancers_1248_ACTIVE.simple.NORM$minus)
median(enhancers_1248_ACTIVE.simple.NORM$plus)
median(enhancers_1248_ACTIVE.simple.NORM$minus1)
median(enhancers_1248_ACTIVE.simple.NORM$plus1)
     
median(enhancers_LESS_ACTIVE.simple.NORM$minus) 
median(enhancers_LESS_ACTIVE.simple.NORM$plus)
median(enhancers_LESS_ACTIVE.simple.NORM$minus1)
median(enhancers_LESS_ACTIVE.simple.NORM$plus1)

############################################################################################

wilcox.test(enhancers_1248_ACTIVE.simple.NORM$minus,   enhancers_1248_ACTIVE.simple.NORM$plus)$p.value
wilcox.test(enhancers_1248_ACTIVE.simple.NORM$minus1,  enhancers_1248_ACTIVE.simple.NORM$plus1)$p.value 
wilcox.test(enhancers_1248_ACTIVE.simple.NORM$plus,    enhancers_1248_ACTIVE.simple.NORM$plus1) 

############################################################################################

wilcox.test(enhancers_LESS_ACTIVE.simple.NORM$minus, enhancers_LESS_ACTIVE.simple.NORM$plus)$p.value  
wilcox.test(enhancers_LESS_ACTIVE.simple.NORM$minus1, enhancers_LESS_ACTIVE.simple.NORM$plus1)$p.value  
wilcox.test(enhancers_LESS_ACTIVE.simple.NORM$plus,  enhancers_LESS_ACTIVE.simple.NORM$plus1)$p.value   

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ the BOXPLOTS of 4 

TITLE_HERE="comp:enhancers:1248:OTHERS.enhancers"
COLORS2 = c("red", "blue", "red", "blue")
LIMITS=c(0,400)

png(paste(TITLE_HERE, "display.the.BOXPLOTS.of.CPM", "png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(enhancers_1248_ACTIVE.simple.NORM$minus, 
        enhancers_1248_ACTIVE.simple.NORM$plus,
        enhancers_1248_ACTIVE.simple.NORM$minus1,
        enhancers_1248_ACTIVE.simple.NORM$plus1,

        enhancers_OTHERS.simple.NORM$minus, 
        enhancers_OTHERS.simple.NORM$plus,
        enhancers_OTHERS.simple.NORM$minus1,
        enhancers_OTHERS.simple.NORM$plus1,

        ylim=LIMITS, 
        col=COLORS2, 
        notch=TRUE, 
        main=TITLE_HERE, 
        # xlab=c("minus","plus"),   
        ylab="normalized counts")
dev.off()  

############################################################################################

median(enhancers_1248_ACTIVE.simple.NORM$minus)
median(enhancers_1248_ACTIVE.simple.NORM$plus)
median(enhancers_1248_ACTIVE.simple.NORM$minus1)
median(enhancers_1248_ACTIVE.simple.NORM$plus1)

median(enhancers_OTHERS.simple.NORM$minus)
median(enhancers_OTHERS.simple.NORM$plus)
median(enhancers_OTHERS.simple.NORM$minus1)
median(enhancers_OTHERS.simple.NORM$plus1)

############################################################################################

wilcox.test(enhancers_1248_ACTIVE.simple.NORM$minus, enhancers_1248_ACTIVE.simple.NORM$plus)$p.value 
wilcox.test(enhancers_1248_ACTIVE.simple.NORM$minus1, enhancers_1248_ACTIVE.simple.NORM$plus1)$p.value 
wilcox.test(enhancers_1248_ACTIVE.simple.NORM$plus, enhancers_1248_ACTIVE.simple.NORM$plus1)$p.value 

wilcox.test(enhancers_OTHERS.simple.NORM$minus, enhancers_OTHERS.simple.NORM$plus)$p.value 
wilcox.test(enhancers_OTHERS.simple.NORM$minus1, enhancers_OTHERS.simple.NORM$plus1)$p.value 
wilcox.test(enhancers_OTHERS.simple.NORM$plus, enhancers_OTHERS.simple.NORM$plus1)$p.value 

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

TITLE_HERE="comp:here:the.MEDIANS:and:the.PVALUES.txt"

############################################################################################

write.table("enhancers_1248_ACTIVE.simple.NORM$minus, enhancers_1248_ACTIVE.simple.NORM$plus \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table(wilcox.test(enhancers_1248_ACTIVE.simple.NORM$minus, enhancers_1248_ACTIVE.simple.NORM$plus)$p.value,
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table("enhancers_1248_ACTIVE.simple.NORM$minus1, enhancers_1248_ACTIVE.simple.NORM$plus1 \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table(wilcox.test(enhancers_1248_ACTIVE.simple.NORM$minus1, enhancers_1248_ACTIVE.simple.NORM$plus1)$p.value,
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table("enhancers_1248_ACTIVE.simple.NORM$plus, enhancers_1248_ACTIVE.simple.NORM$plus1 \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table(wilcox.test(enhancers_1248_ACTIVE.simple.NORM$plus, enhancers_1248_ACTIVE.simple.NORM$plus1)$p.value,
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

############################################################################################

write.table("enhancers_LESS_ACTIVE.simple.NORM$minus, enhancers_LESS_ACTIVE.simple.NORM$plus \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table(wilcox.test(enhancers_LESS_ACTIVE.simple.NORM$minus, enhancers_LESS_ACTIVE.simple.NORM$plus)$p.value,
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table("enhancers_LESS_ACTIVE.simple.NORM$minus1, enhancers_LESS_ACTIVE.simple.NORM$plus1 \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table(wilcox.test(enhancers_LESS_ACTIVE.simple.NORM$minus1, enhancers_LESS_ACTIVE.simple.NORM$plus1)$p.value,
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table("enhancers_LESS_ACTIVE.simple.NORM$plus, enhancers_LESS_ACTIVE.simple.NORM$plus1 \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table(wilcox.test(enhancers_LESS_ACTIVE.simple.NORM$plus, enhancers_LESS_ACTIVE.simple.NORM$plus1)$p.value,
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

############################################################################################ 

write.table("enhancers_OTHERS.simple.NORM$minus, enhancers_OTHERS.simple.NORM$plus \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table(wilcox.test(enhancers_OTHERS.simple.NORM$minus, enhancers_OTHERS.simple.NORM$plus)$p.value,
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 


write.table("enhancers_OTHERS.simple.NORM$minus1, enhancers_OTHERS.simple.NORM$plus1 \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(wilcox.test(enhancers_OTHERS.simple.NORM$minus1, enhancers_OTHERS.simple.NORM$plus1)$p.value,
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table("enhancers_OTHERS.simple.NORM$plus, enhancers_OTHERS.simple.NORM$plus1 \n",
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

write.table(wilcox.test(enhancers_OTHERS.simple.NORM$plus, enhancers_OTHERS.simple.NORM$plus1)$p.value, 
            file=TITLE_HERE, 
            append=TRUE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE) 

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################ 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
