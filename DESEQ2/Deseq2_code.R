#############DESEQ2 for info.csv##############################

setwd("/Users/achisha_saikia/Documents/Salivary glands/Mouse salivary glands/DESEQ2_All_glands_CD1_vs_C57/")
library("DESeq2")
library(ggplot2)
library(ggrepel)
library("RColorBrewer")
library("pheatmap")
library(cowplot)
library('genefilter')
library('gplots')
library("limma")
library(apeglm)
library(EnhancedVolcano)
library(dplyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(genefilter)
library(PoiClaClu)

#Create the dds object with the original featurecounts and include batch in your design formula. Your design formula may look like ~ batch + condition 
#Conduct a differential expression analysis in the normal way for your factor of interest (e.g., condition) - the effect of batch will be 'adjusted'
#when the statistical inferences are made due to the fact that batch is included in the design formula
#Run the commands with limma::removeBatchEffect() if you want to use your data matrix downstream for clustering, PCA, machine learning
#########################################################################################################################################################
# Read in counts table 'cts' from featureCount output
getwd()
setwd("/Users/achisha_saikia/Documents/Salivary glands/Mouse salivary glands/DESEQ2_All_glands_CD1_vs_C57/")
dat <- read.csv("readCounts.tsv",
                header = TRUE, 
                sep = '\t',
                row.names = 'Geneid',
                skip=1)

# Drop extra columns from featureCount output
drops <- c("Chr","Start","End","Strand","Length")
head(drops)
dat <- dat[ , !(names(dat) %in% drops)]
colnames(dat)


#If I want to subset instead of dropping columns
#colnames(dat) <- c("Par_1A_CD1","Par_1B_CD57","Par_2A_CD1","Par_2B_CD57","Par_3A_CD1","Par_3B_CD57")

# Read in sample information table 'coldata'
info <- read.csv("info.csv", header = TRUE, 
                 sep = ',',
                 stringsAsFactors = TRUE)

dds <- DESeqDataSetFromMatrix(countData = dat, colData = info, design = ~ Condition)

#remove lowly expressed genes
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]

#main DESeq
dds<-DESeq(dds)

#Compare and contrast
result<- results(dds, contrast = c("Condition", "C57", "CD1"))
resultsNames(dds)

#res1_B6_vs_CD1 <- results(dds, contrast = list("Condition_B6_vs_CD1"))
#res1_B6_vs_CD57 <- results(dds, contrast = list("Condition_B6_vs_CD57"))
#res1_CD1_vs_CD57 <- results(dds, contrast = list("Condition_B6_vs_CD1", "Condition_B6_vs_CD57"))
#resultsNames(ddsDE)

#export normalized read counts
normCounts <- counts(dds, normalized=T)
write.csv(normCounts, "normalized_All_glands_CD1_vs_C57_counts.csv")

#DESeq results
res<-results(dds, alpha = 0.05)

resOrdered <- res[order(res$padj),]
write.csv(resOrdered, "deseq_All_glands_CD1_vs_C57_counts.csv")

summary(res)
resultsNames(dds)



#####If I want to load the files again######################
normCount <- read.csv("normalized_All_glands_CD1_vs_C57_counts.csv", row.names = 1)
deSeqRes <- read.csv("deseq_All_glands_CD1_vs_C57_counts.csv", row.names = 1)
View(deSeqRes)


#Setting parameters
#deSeqRes$sig <- ifelse((deSeqRes$log2FoldChange > 1.5 & deSeqRes$padj <0.05), "yes", "no")
deSeqRes$sig <- ifelse(deSeqRes$padj <= 0.05, "yes", "no")
deSeqRes <- na.omit(deSeqRes)

write.csv(deSeqRes, "Deseq2_significant_res_All_glands_CD1_vs_C57.csv")

View(deSeqRes)


##########Plotting##################

#MA Plot
ggplot(deSeqRes, aes(x= log10(baseMean), y= log2FoldChange, color = sig)) + geom_point()


#volcana plot
ggplot(deSeqRes, aes(x= log2FoldChange, y = -log10(padj), color= sig)) + geom_point() + scale_color_manual(values=c("maroon","blue")) + theme_minimal() 


#pheatmap
signi <- subset(deSeqRes, padj <= 0.05)
#View(signi)
allsig <- merge(normCount, signi, by = 0)
#View(allsig)

sigCounts <- allsig[,2:36]
row.names(sigCounts) <- allsig$Row.names
#View(sigCounts)
pheatmap(log2(sigCounts +1), scale = 'row', show_rownames = F, treeheight_row = 0, treeheight_col = 0)


#PCA
vsdata<-vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "Sample") 
plotDispEsts(dds)


######counts data

d <- plotCounts(dds, gene="Pgam2", intgroup=c("Sample","Condition", "Design"), 
                returnData=TRUE)
ggplot(d, aes(x=Sample, y=count,color=Condition,shape=Condition)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + theme_minimal() + ggtitle("Pgam2")







