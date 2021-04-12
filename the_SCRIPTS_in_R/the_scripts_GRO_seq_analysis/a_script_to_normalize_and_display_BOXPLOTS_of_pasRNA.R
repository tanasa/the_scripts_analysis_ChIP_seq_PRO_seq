##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

library("edgeR")
library("ggplot2")
library("extrafont")
loadfonts()

##############################################################################################################
##############################################################################################################
##############################################################################################################
#################################################################### reading the files with GENE ID, and NAMES

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

genes.RAND = read.delim("genes.random837.bed",
                       header = FALSE, sep = "\t", stringsAsFactors=F)
dim(genes.RAND)
head(genes.RAND)

####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

genes.all.simple   = data.frame(RefSeq=genes.all$V4)
genes.PR.simple    = data.frame(RefSeq=genes.PR$V4)
genes.OG.simple    = data.frame(RefSeq=genes.OG$V4)
genes.RAND.simple  = data.frame(RefSeq=genes.RAND$V4)

dim(genes.all.simple) 
dim(genes.PR.simple)  
dim(genes.OG.simple)  
dim(genes.RAND.simple)  

##############################################################################################################
##############################################################################################################
##############################################################################################################
#################################################################### here reading the file with GRO-seq
#################################################################### here the values of RAW COUNTS

FILE = ""  ############ reading the file with the GRO-seq RAW COUNTS 

x <- read.delim(FILE,
                 header=T,
                 sep="\t", stringsAsFactors=F)

dim(x)
head(x)

x_name <- "pasRNA.hg18.counts.RAW.size.1kb"
CHIP="pasRNA"
colnames(x)

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
##################################################################################### we can select SPECIFIC COLUMNS 

x.small <- data.frame(RefSeq = x$REFSEQ,
                      gene = x$GENE,
                      minus = x$OPPOSITE_STRAND_0m_TSS,
                      plus = x$OPPOSITE_STRAND_40m_TSS)    

dim(x.small)  
head(x.small)

#####################################################################################################################
################################################################################### to work with the LONGEST ISOFORMS

dim(genes.all.simple) 
dim(genes.PR.simple)  
dim(genes.OG.simple)  
dim(genes.RAND.simple) 

x.small.LI <- merge(genes.all.simple,
                    x.small,
                    by.x = "RefSeq",
                    by.y = "gene",
                    all = FALSE )

dim(x.small.LI) 
head(x.small.LI)

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################## a dataframe for CPM and NORMALIZATION: edgeR 

row.names(x.small.LI) <- x.small.LI$RefSeq

x.small.LI <- x.small.LI[,-1]
head(x.small.LI)
x.small.LI <- x.small.LI[,-1]
head(x.small.LI)

#######################################################################################################################
#######################################################################################################################

group <- factor(c("minus", "plus"))
libSizes <- as.vector(colSums(x.small.LI))
y <- DGEList(counts=x.small.LI, group=group, lib.size=libSizes)

#######################################################################################################################
#######################################################################################################################
####################################################################################################################### 
#######################################################################################################################

# y <- calcNormFactors(y)
# counts.per.m <- cpm(y, normalized.lib.sizes=TRUE)

#######################################################################################################################
############################################################################################################## CPM > 1
#######################################################################################################################

z <- y  ##### z to be used when doing the filtering 
keep <- rowSums(cpm(z) > 1) >= 1
z <- z[keep, ,keep.lib.sizes=FALSE]
dim(z)

#######################################################################################################################
####################################################################################################### NORMALIZATION :
#######################################################################################################################

z <- calcNormFactors(z)

counts.per.m.z <- cpm(z, normalized.lib.sizes=TRUE)
counts.per.m.z <- as.data.frame(counts.per.m.z)

png(paste(x_name, "ALL_genes", "display_after_NORMALIZATION.edgeR.filtering.png", sep="."))
boxplot(counts.per.m.z[,1], 
        counts.per.m.z[,2],
        ylim=c(0,200),
        col="yellow", main=" NORMALIZED COUNTS in edgeR after FILTERING")
dev.off()

# wilcox.test(counts.per.m.z$minus, counts.per.m.z$plus)

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
############################################################################################################## here to start using : counts.per.m.z 
############################################################################################################## and to integrate with the DATAFRAMES

counts.per.m.z <- as.data.frame(counts.per.m.z)
counts.per.m.z$gene <- row.names(counts.per.m.z) 

head(counts.per.m.z)
dim(counts.per.m.z)

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
############################################################################################################## using counts.per.m.z

length(genes.PR.simple$RefSeq)  
length(genes.OG.simple$RefSeq)  
length(genes.RAND.simple$RefSeq) 

##############################################################################################################
##############################################################################################################
##############################################################################################################
############################################################################################################## genes.PR.simple

genes.PR.simple.NORM <- merge(genes.PR.simple,
                              counts.per.m.z,
                              by.x = "RefSeq",
                              by.y = "gene",
                              all = FALSE )

dim(genes.PR.simple) 
dim(genes.PR.simple.NORM)  
dim(counts.per.m.z)    
head(genes.PR.simple.NORM)

#######################################################################
#######################################################################

png(paste(x_name, "counts.per.m.z.genes.PR.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.NORM$minus,
        genes.PR.simple.NORM$plus,
        ylim=c(0,200),
        col=c("red", "blue"), 
        notch=TRUE, 
        main="genes.UP.PR.edgeR.norm")
dev.off()

# wilcox.test(genes.PR.simple.NORM$minus, genes.PR.simple.NORM$plus)

##############################################################################################################
##############################################################################################################
##############################################################################################################
############################################################################################################## genes.OG.simple

genes.OG.simple.NORM <- merge(genes.OG.simple,
                              counts.per.m.z,
                              by.x = "RefSeq",
                              by.y = "gene",
                              all = FALSE )

dim(genes.OG.simple)      
dim(genes.OG.simple.NORM)  
dim(counts.per.m.z)       
head(genes.OG.simple.NORM)

#######################################################################
#######################################################################

png(paste(x_name, "counts.per.m.z.genes.OG.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.OG.simple.NORM$minus,
        genes.OG.simple.NORM$plus,
        ylim=c(0,200),
        col=c("red", "blue"), 
        notch=TRUE, 
        main="genes.OG.edgeR.norm")
dev.off()

# wilcox.test(genes.OG.simple.NORM$minus,genes.OG.simple.NORM$plus)

##############################################################################################################
##############################################################################################################
##############################################################################################################
############################################################################################################## genes.RAND.simple 

genes.RAND.simple.NORM <- merge(genes.RAND.simple,
                                counts.per.m.z,
                                by.x = "RefSeq",
                                by.y = "gene",
                                all = FALSE )

dim(genes.RAND.simple)      
dim(genes.RAND.simple.NORM) 
dim(counts.per.m.z)         
head(genes.RAND.simple.NORM)

#######################################################################
#######################################################################

png(paste(x_name, "counts.per.m.z.genes.RAND.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.RAND.simple.NORM$minus,
        genes.RAND.simple.NORM$plus,
        ylim=c(0,200),
        col=c("red", "blue"), 
        notch=TRUE, 
        main="genes.RAND.edgeR.norm")
dev.off()

# wilcox.test(genes.RAND.simple.NORM$minus, genes.RAND.simple.NORM$plus)

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
###################################################### here making the BOXPLOTS for each of these CATEGORIES :

#genes.PR.simple.NORM$minus,
#genes.PR.simple.NORM$plus,

#genes.OG.simple.NORM$minus,
#genes.OG.simple.NORM$plus,

#genes.RAND.simple.NORM$minus,
#genes.RAND.simple.NORM$plus,

#counts.per.m.z$minus,
#counts.per.m.z$plus,

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

COLORS=c("green", "green", "red", "red", "magenta", "magenta", "cyan", "cyan", "blue", "blue")

COLORS_for_OG=c("red", "blue", "red", "blue")

COLORS_for_RAND=c("red", "blue", "red", "blue")

##############################################################################################################
############################################################################################################## with OG

png(paste(x_name, "counts.per.m.genes.2.classes.UP_PR.UP_noPR.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.NORM$minus,
        genes.PR.simple.NORM$plus,
        genes.OG.simple.NORM$minus,
        genes.OG.simple.NORM$plus,
        #genes.RAND.simple.NORM$minus,
        #genes.RAND.simple.NORM$plus,
        #counts.per.m.z$minus,
        #counts.per.m.z$plus,
        ylim=c(0,200),  
        col=COLORS_for_OG, 
        notch=TRUE, 
        main="UP_PR::UP_noPR(OG)", 
        ylab="normalized counts")
dev.off()

##############################################################################################################
############################################################################################################## with RAND

png(paste(x_name, "counts.per.m.genes.2.classes.UP_PR.RAND.png", sep="."), height=12, width=8, units="cm", res=800)
boxplot(genes.PR.simple.NORM$minus,
        genes.PR.simple.NORM$plus,
        #genes.OG.simple.NORM$minus,
        #genes.OG.simple.NORM$plus,
        genes.RAND.simple.NORM$minus,
        genes.RAND.simple.NORM$plus,
        #counts.per.m.z$minus,
        #counts.per.m.z$plus,
        ylim=c(0,200),  
        col=COLORS_for_RAND, 
        notch=TRUE, 
        main="UP_PR::RAND", 
        ylab="normalized counts")
dev.off()

#########################################################################################
#########################################################################################
######################################################################################### to OUTPUT in a FILE
#########################################################################################
#########################################################################################

median(genes.PR.simple.NORM$minus)
median(genes.PR.simple.NORM$plus)
median(genes.OG.simple.NORM$minus)
median(genes.OG.simple.NORM$plus)
median(genes.RAND.simple.NORM$minus)
median(genes.RAND.simple.NORM$plus)
median(counts.per.m.z$minus)
median(counts.per.m.z$plus)

######################################################################################### 
#########################################################################################
######################################################################################### the P-values
######################################################################################### wilcoxon.test

wilcox.test(counts.per.m.z$minus, counts.per.m.z$plus)$p.value   

wilcox.test(genes.PR.simple.NORM$minus, genes.PR.simple.NORM$plus)$p.value  

wilcox.test(genes.OG.simple.NORM$minus, genes.OG.simple.NORM$plus)$p.value  

wilcox.test(genes.RAND.simple.NORM$minus, genes.RAND.simple.NORM$plus)$p.value 

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

write.table("the RESULTS of WILCOXON TESTS : ALL genes", 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(wilcox.test(counts.per.m.z$minus, counts.per.m.z$plus)$p.value, 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table("the RESULTS of WILCOXON TESTS : UP_PR", 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(wilcox.test(genes.PR.simple.NORM$minus, genes.PR.simple.NORM$plus)$p.value, 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table("the RESULTS of WILCOXON TESTS : UPnonPR", 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(wilcox.test(genes.OG.simple.NORM$minus, genes.OG.simple.NORM$plus)$p.value, 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table("the RESULTS of WILCOXON TESTS : RAND2", 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(wilcox.test(genes.RAND.simple.NORM$minus, genes.RAND.simple.NORM$plus)$p.value, 
            file = "the.results.WILCOX.test.txt", 
            append = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
