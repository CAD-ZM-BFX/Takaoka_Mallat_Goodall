## /local/usr/bin/Rscript
## CAD_jcg23_0001 
## RNASeq data mouse, 4 time points (0, 2, 4, 7), KO vs WT.
## 11/2021
# Analysis Performed by Xiaohui Zhao
# School of Medicine, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+---------------- Basic libraries and settings loading -------------------+")
suppressPackageStartupMessages({
  library('DESeq2')
  library('ggplot2')
  library('RColorBrewer')
  library("cowplot")
  library("reshape")
  library("pheatmap")
  library("ggrepel")
  library("reshape2")
  library("biomaRt")
  library("matrixStats")
  library("plyr")
  library("BiocParallel")
  library("dplyr")
  library("ggalt")
  library("limma")
  library("apeglm")
  library("gdata") 
  library("ComplexHeatmap")
  library("seriation")
  library("methods")
  library("utils")
  library("Matrix")
  library("useful")
  library("edgeR")
  library("GetoptLong")
  library("UpSetR")
  library("circlize")
  library("VennDiagram")
  library("eulerr")
  library("readr")
  library("gridExtra")
  library("grid")
  library("hexbin")
  library("xlsx")
  ## GO package
  library("clusterProfiler")
  library("DOSE")
  library("GSEABase")
  library("AnnotationHub")
  library("org.Mm.eg.db")
  library("gage")
  library("gageData")
  library("enrichplot")
  library("ggraph")
  library("ggforce")
  library("vsn")
  ## secreted and tissueEnrich
  library("TissueEnrich")
  library("UniProt.ws")
  library("readxl")
  library("corrplot") 
  library("clusterSim") 
})
register(MulticoreParam(2))
setwd("/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001")
Project         <- "CAD_jcg23_0001"
significance    <- 0.05
l2fc            <- 1 
elementTextSize <- 10
TOPNUM          <- 2000
elementTextSize <- 10

## load ensembl annotation GRCm39_104
ensembl    =  useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host = 'ensembl.org')
listEnsembl()
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description','chromosome_name',
                                 'start_position', 'end_position', 'strand', 'transcript_length'),  mart = ensembl, useCache = FALSE) 
ensEMBL2id$description <- gsub("..Source.*", "", ensEMBL2id$description)
ensEMBL2id <- ensEMBL2id[-which(duplicated(ensEMBL2id$ensembl_gene_id)==T),] 
## Keep the longest transcipt length.
ensEMBL2id <- subset(ensEMBL2id, chromosome_name==1|chromosome_name==2|chromosome_name==3|chromosome_name==4|
                       chromosome_name==5|chromosome_name==6|chromosome_name==7|chromosome_name==8|
                       chromosome_name==9|chromosome_name==10|chromosome_name==11|chromosome_name==12|
                       chromosome_name==13|chromosome_name==14|chromosome_name==15|chromosome_name==16|
                       chromosome_name==17|chromosome_name==18|chromosome_name==19|chromosome_name=="MT"|
                       chromosome_name=="X"|chromosome_name=="Y")

head(ensEMBL2id)
nrow(ensEMBL2id) ## 55359

save(ensEMBL2id, file = "./Original_Data/Ensembl_mmusculus_ID_Name_Des_Chr_GRCm39_ensembl_104.RData")
#load("./Original_Data/Ensembl_mmusculus_ID_Name_Des_Chr_GRCm39_ensembl_104.RData")

message("+---- Calling the gene sample count matrix and sample Table --------------------------------+")
 
count_df_all       <- read.table("./Original_Data/CAD_jcg23_0001-RawCounts.txt", header=T)
colnames(count_df_all)[4:6] <- c("WT_T2w_REP4", "KO_T2w_REP1", "KO_T2w_REP2")
## Jane figured out that in 2w one KO is wrong with label, it should be the WT,
count_df_all       <- count_df_all[,order(colnames(count_df_all))]
targets_all        <- as.data.frame(read_excel("./Original_Data/SLX-21013_Summary_Table.xlsx", sheet = 1))
colnames(targets_all) <- c("SampleName", "Condition", "Time", "Sex", "GDF15")
targets_all$Time   <- ifelse(targets_all$Time=="Non-irradiated", "0w", targets_all$Time)
targets_all$Time   <- ifelse(targets_all$Time=="2 weeks", "2w", targets_all$Time)
targets_all$Time   <- ifelse(targets_all$Time=="4 weeks", "4w", targets_all$Time)
targets_all$Time   <- ifelse(targets_all$Time=="7 weeks", "7w", targets_all$Time)
targets_all$Rep    <- c(rep(c(1:3), length = 6), c(1:4), c(1:2), rep(c(1:3), length=12))
targets_all$sample <- paste0(targets_all$Condition, "_T", targets_all$Time, "_REP", targets_all$Rep)
rownames(targets_all) <- targets_all$sample
targets_all        <- targets_all[order(targets_all$sample),]

message("+---------------- Perform DESeq analysis as interested            --------------------------------+")

dds_all   <- DESeqDataSetFromMatrix(countData=count_df_all, colData=targets_all, design=~Time+Condition)
dds_all$Condition <- relevel(dds_all$Condition, ref = "WT")
dds_all   <- DESeq(dds_all, parallel=TRUE)
dds_all   <- estimateSizeFactors(dds_all)
vsd_all   <- vst(dds_all,     blind=F)
colData(vsd_all)
resultsNames(dds_all)
normalized_counts <- as.data.frame(counts(dds_all, normalized=TRUE))
normalized_counts$ensembl_gene_id <- rownames(normalized_counts)
normcounts <- merge(normalized_counts, ensEMBL2id, by ="ensembl_gene_id")
write.csv(normcounts, file = "./Results_Tables_Figures/March_2023/CAD_jcg23_0001-AllSamples_Normalised_Counts_March_2023.csv")

message("+--------------- Informal PCA plot function and analysis -----------------------+")
customPCA  <- function(sampleTBL, RLD, TOPNUM, model, ensEMBL2id) {
  RLD                 <- as.data.frame(RLD) 
  RLD.ori             <- RLD
  RLD.ori$ensembl_gene_id <- rownames(RLD.ori)
  RLD.mer             <- merge(RLD.ori, ensEMBL2id, by="ensembl_gene_id") 
  RLD.mer             <- RLD.mer[-which(RLD.mer$external_gene_name==""),]## remove no gene names
  RLD.mer             <- RLD.mer[-(which(duplicated(RLD.mer$external_gene_name)==T)),] # remove duplicated genes
  RLD.new             <- RLD.mer[,colnames(RLD)]
  rownames(RLD.new)   <- RLD.mer$external_gene_name
  colnames(RLD.new)   <- sampleTBL$SampleName
  rv     <- rowVars(as.matrix(RLD.new))
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD.new[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  rownames(sampleTBL) <- sampleTBL$SampleName
  scores.dat <- merge(sampleTBL, pca$x, by="row.names")
  
  scores    <- data.frame(sampleName=scores.dat$sample, scores.dat[,c("PC1", "PC2")], 
                          condition=scores.dat$Condition,
                          time=gsub("w", "", scores.dat$Time),
                          condition_week=paste0(scores.dat$Condition,"_", gsub("w", "", scores.dat$Time)))
  
  scores.summary           <- ddply(scores, c("condition_week"), summarise, meanPC1 = median(PC1), meanPC2 = median(PC2))
  scores.summary$week      <- scores.summary$condition_week
  scores.summary$week      <- gsub(".*_", "",scores.summary$week)
  scores.summary$week      <- as.integer(scores.summary$week)
  
  scores.summary           <- scores.summary[order(scores.summary$week, decreasing = FALSE),]
  scores.summary$condition <- scores.summary$condition_week
  scores.summary$condition <- substr(scores.summary$condition, 1,2)
  scores.summary$condition <- factor(scores.summary$condition, levels = c("WT", "KO"))
  scores.summary$week      <- factor(scores.summary$week, levels = c("0","2", "4", "7"))
  
  
  
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=condition_week) ) +
    geom_point(size = 3, alpha=0.75) + 
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
    theme(text = element_text(size=elementTextSize)) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, label=condition, group=condition_week) ) +
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(fill=condition)) +
    geom_point(size = 3, alpha=0.75 ) +
    geom_text_repel(data=scores.summary, aes(x=meanPC1, y=meanPC2, label=condition_week), show.legend = FALSE, size=4, colour="black") +
    scale_colour_manual(name="Treatment", values=c("WT"="royalblue",  "KO"="darkorange")) +
    scale_fill_manual(name="Treatment", values=c("WT"="royalblue", "KO"="darkorange")) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize), legend.position="right") +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  plt.pca.nl_path <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, label=condition, group=condition_week) ) +
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(fill=condition)) +
    
    geom_point(data=scores.summary, aes(x=meanPC1, y=meanPC2, group=condition, colour=condition), shape=5, size=5, alpha=0.75, show.legend=TRUE) + 
    geom_path(data=scores.summary,  aes(x=meanPC1, y=meanPC2, group=condition), show.legend=FALSE)  + 
    
    geom_point(size = 3, alpha=0.75 ) + 
    geom_text_repel(data=scores.summary, aes(x=meanPC1, y=meanPC2, label=condition_week), show.legend = FALSE, size=4, colour="black") +
    
    
    scale_colour_manual(name="Treatment", values=c("WT"="royalblue",  "KO"="darkorange")) +
    scale_fill_manual(name="Treatment", values=c("WT"="royalblue", "KO"="darkorange")) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize), legend.position="right") +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  loadings                 <- as.data.frame(pca$rotation)
  
  
  pca.1         <- loadings[order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:25),]
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(rownames(pca.1.25),levels=unique(rownames(pca.1.25))), y=PC1)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  pca.2         <- loadings[order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:25),]
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(rownames(pca.2.25),levels=unique(rownames(pca.2.25))), y=PC2)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme()+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  return(list(plt.pca,plt.pca.nl,plt.pca.nl_path,pca.1.25.plot, pca.2.25.plot))
  
}

pca_all    <- customPCA(targets_all, assay(vsd_all), TOPNUM, "GADD34", ensEMBL2id)

pdf("./Results_Tables_Figures/March_2023/CAD_jcg23_0001-Cust_PCA_allSamples_March_2023.pdf")
print(pca_all[[2]])
dev.off()
## modified adding the path following time and also change the coloron 24/04/2023
pdf("./Results_Tables_Figures/Paper_v1_April_2023/SuppFigxx-Cust_PCA_allSamples_April_2023_withoutPath.pdf")
print(pca_all[[2]])
dev.off()
pdf("./Results_Tables_Figures/Paper_v1_April_2023/SuppFigxx-Cust_PCA_allSamples_April_2023_withPath.pdf")
print(pca_all[[3]])
dev.off()

message("+----------- WT DESeq analysis across time points-------------------------+")

counts_df_WT <- count_df_all[,grep("WT", colnames(count_df_all))]
targets_WT   <- targets_all[grep("WT", targets_all$Condition),]
dds_WT  <- DESeqDataSetFromMatrix(countData=counts_df_WT, colData=targets_WT, design=~Time)
#dds_WT$Time <- relevel(dds_WT$Time, ref = "0w")
dds_WT   <- DESeq(dds_WT, parallel=TRUE)
dds_WT   <- estimateSizeFactors(dds_WT)
vsd_WT   <- vst(dds_WT,     blind=F)
colData(vsd_WT)
resultsNames(dds_WT)

pca_WT    <- customPCA(targets_WT, assay(vsd_WT), TOPNUM, "GADD34_WT", ensEMBL2id)

pdf("./Results_Tables_Figures/March_2023/WT_allsamples_T0247_PCAplot_March_2023.pdf")
plot_grid(pca_WT[[1]], pca_WT[[2]], pca_WT[[3]], pca_WT[[4]], nrow=2, byrow=T)
dev.off()

message("+----------- KO DESeq analysis across time points-------------------------+")

counts_df_KO <- count_df_all[,grep("KO", colnames(count_df_all))]
targets_KO   <- targets_all[grep("KO", targets_all$Condition),]
dds_KO  <- DESeqDataSetFromMatrix(countData=counts_df_KO, colData=targets_KO, design=~Time)
#dds_WT$Time <- relevel(dds_WT$Time, ref = "0w")
dds_KO   <- DESeq(dds_KO, parallel=TRUE)
dds_KO   <- estimateSizeFactors(dds_KO)
vsd_KO   <- vst(dds_KO,     blind=F)
colData(vsd_KO)
resultsNames(dds_KO)

pca_KO    <- customPCA(targets_KO, assay(vsd_KO), TOPNUM, "GADD34_KO", ensEMBL2id)

pdf("./Results_Tables_Figures/March_2023/KO_allsamples_T0247_PCAplot_March_2023.pdf")
plot_grid(pca_KO[[1]], pca_KO[[2]], pca_KO[[3]], pca_KO[[4]], nrow=2, byrow=T)
dev.off()

message("+--- Formal hierchary clustering analysis plot for the all samples and WT samples -------------+")

rv_all  <- rowVars(assay(vsd_all))
o_all   <- order(rv_all,decreasing=TRUE)
dists_all <- dist(t(assay(vsd_all)[head(o_all,TOPNUM),]))
hc_all    <- hclust(dists_all)
pdf("./Results_Tables_Figures/March_2023/Hclust_March_2023.pdf", width=8, height=5)
plot(hc_all, labels=vsd_all$sample)
dev.off()

rv_WT  <- rowVars(assay(vsd_WT))
o_WT   <- order(rv_WT,decreasing=TRUE)
dists_WT <- dist(t(assay(vsd_WT)[head(o_WT,TOPNUM),]))
hc_WT    <- hclust(dists_WT)
pdf("./Results_Tables_Figures/March_2023/Hclust_WTonly_March_2023.pdf", width=8, height=5)
plot(hc_WT, labels=vsd_WT$sample)
dev.off()

rv_KO  <- rowVars(assay(vsd_KO))
o_KO   <- order(rv_KO,decreasing=TRUE)
dists_KO <- dist(t(assay(vsd_KO)[head(o_KO,TOPNUM),]))
hc_KO    <- hclust(dists_KO)
pdf("./Results_Tables_Figures/March_2023/Hclust_KOonly_March_2023.pdf", width=8, height=5)
plot(hc_KO, labels=vsd_KO$sample)
dev.off()

message("+----Paired DEGs analysis for each time points. ---------------------+")

customPCA_pair     <- function(sampleTBL, ensEMBL2id, RLD, TOPNUM) {
  RLD            <- as.data.frame(RLD) 
  RLD$ensembl_gene_id <- rownames(RLD)
  RLD.mer<- merge(RLD, ensEMBL2id, by = "ensembl_gene_id")
  RLD.mer<- RLD.mer[-which(duplicated(RLD.mer$external_gene_name)==T),]
  rownames(RLD.mer) <- RLD.mer$external_gene_name
  RLD.mer<- RLD.mer[,2:7]
  rv     <- rowVars(as.matrix(RLD.mer))
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD.mer[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sample, pca$x, 
                          condition=sampleTBL$Condition)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3, alpha=0.75) + 
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
    theme(text = element_text(size=elementTextSize)) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3, alpha=0.75 ) + 
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(colour=condition, fill=condition, group = condition)) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize))+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  loadings                 <- as.data.frame(pca$rotation)
  pca.1         <- loadings[order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:25),]
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(rownames(pca.1.25),levels=unique(rownames(pca.1.25))), y=PC1)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  pca.2         <- loadings[order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:25),]
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(rownames(pca.2.25),levels=unique(rownames(pca.2.25))), y=PC2)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme()+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  return(list(plt.pca,plt.pca.nl,pca.1.25.plot, pca.2.25.plot))
  
}
merECRes_function_pair <- function(resfile, ensEMBL2id, ddsfile, outfile){
  resdat   <- as.data.frame(resfile)
  resdat$ensembl_gene_id <- rownames(resdat)
  resdat.merE <- merge(resdat, ensEMBL2id, by ="ensembl_gene_id")
  normCnt     <- as.data.frame(counts(ddsfile, normalized=TRUE))
  normCnt$ensembl_gene_id <- rownames(normCnt)
  resmerEC <- merge(resdat.merE, normCnt, by = "ensembl_gene_id")
  resmerEC <- subset(resmerEC, !is.na(padj) & !is.na(log2FoldChange))
  print(dim(resmerEC)) 
  write.csv(resmerEC, file=outfile, row.names=F)
}

Time <- c("0w", "2w", "4w", "7w")

RESsig.files <- paste0("./Results_Tables_Figures/March_2023/",Project, "-ResSig_",Time, "_summaryTable_March_2023.csv")
RES.files <- paste0("./Results_Tables_Figures/March_2023/",Project, "-Resall_",Time, "_summaryTable_March_2023.csv")
RD.files  <- paste0("./Results_Tables_Figures/March_2023/",Project, "-RData_",Time, "_March_2023.RData")
RESsig.files1 <- paste0("./Results_Tables_Figures/March_2023/",Project, "-ResSig_lfcShrink_",Time, "_summaryTable_March_2023.csv")
RES.files1 <- paste0("./Results_Tables_Figures/March_2023/",Project, "-Resall_lfcShrink_",Time, "_summaryTable_March_2023.csv")

for(i in 1:length(Time)){
  count_df <- count_df_all[,grep(Time[i], colnames(count_df_all))]
  targets_sub  <- targets_all[grep(Time[i],targets_all$Time),]
  dds_sub   <- DESeqDataSetFromMatrix(countData=count_df, colData=targets_sub, design=~Condition)
  dds_sub$Condition <- relevel(dds_sub$Condition, ref = "WT")
  dds_sub   <- DESeq(dds_sub, parallel=TRUE)
  dds_sub   <- estimateSizeFactors(dds_sub)
  vsd_sub   <- vst(dds_sub,     blind=F)
  colData(vsd_sub)
  resultsNames(dds_sub)
  normalized_counts <- as.data.frame(counts(dds_sub, normalized=TRUE))
  normalized_counts$ensembl_gene_id <- rownames(normalized_counts)
  normcounts <- merge(normalized_counts, ensEMBL2id, by ="ensembl_gene_id")
  
  pdf(paste0("./Results_Tables_Figures/March_2023/CAD_jcg23_0001-DESeq2_pairwisePCA_", Time[i], "_KOvsWT_March_2023.pdf"), width=8, height= 6)
  pca <- customPCA_pair(targets_sub, ensEMBL2id, assay(vsd_sub), TOPNUM)
  plot_grid(pca[[1]], pca[[2]], pca[[3]], pca[[4]], nrow=2, byrow=T)
  dev.off()
  
  res       <- results(dds_sub, contrast = c("Condition", "KO", "WT"))
  res.sig   <- as.data.frame(subset(res, padj <= significance & abs(log2FoldChange) >= l2fc))
  print(nrow(res.sig)); print(table(sign(res.sig$log2FoldChange)))
  res.sig$ensembl_gene_id <- rownames(res.sig)
  res.sigmer <- merge(res.sig, ensEMBL2id, by ="ensembl_gene_id")
  write.csv(res.sigmer, file = RESsig.files[i])
  res.merEC <- merECRes_function_pair(res, ensEMBL2id, dds_sub, RES.files[i])
  save(dds_sub, vsd_sub, normcounts, file=RD.files[i])
  ## add wald
  res1       <- lfcShrink(dds_sub, coef="Condition_KO_vs_WT", type="apeglm")
  res1.sig   <- as.data.frame(subset(res1, padj <= significance & abs(log2FoldChange) >= l2fc))
  print(nrow(res1.sig)); print(table(sign(res1.sig$log2FoldChange)))
  res1.sig$ensembl_gene_id <- rownames(res1.sig)
  res1.sigmer <- merge(res1.sig, ensEMBL2id, by ="ensembl_gene_id")
  print(nrow(res1.sig)); print(table(sign(res1.sig$log2FoldChange)))
  res1.sig$ensembl_gene_id <- rownames(res1.sig)
  res1.sigmer <- merge(res1.sig, ensEMBL2id, by ="ensembl_gene_id")
  write.csv(res1.sigmer, file = RESsig.files1[i])
  res1.merEC <- merECRes_function_pair(res1, ensEMBL2id, dds_sub, RES.files1[i])
  gc()
  }


message("+----- volcano plot for pairwise analysis between KO vs WT. -----+")

functionPlotDEVolcano <- function(resfile, ensEMBL2id, sig_cut, logfc_cut, title,  xrange, yrange, topN, xlabel, ylabel) {
  
  results       <- read.csv(resfile, header=T)
  results$genes <- results$external_gene_name
  results       <- subset(results, !is.na(padj))
  results       <- results[order(-results$log2FoldChange),]
  
  volc.plt <- ggplot(data=results, aes(x=log2FoldChange, y=-log10(padj), label=genes)) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log10(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    geom_point(data=subset(results, abs(log2FoldChange) < logfc_cut | padj > sig_cut), alpha=0.75, size=0.3, colour="grey") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange >= logfc_cut),      alpha=0.75, size=0.8, colour="red") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange <= -(logfc_cut)),   alpha=0.75, size=0.8, colour="blue") +
    geom_text_repel( data= subset(results, log2FoldChange > logfc_cut & padj<= sig_cut & external_gene_name!="")[1:topN,],
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    geom_text_repel( data= tail(subset(results, log2FoldChange < (-logfc_cut) & padj<= sig_cut & external_gene_name!=""), n=topN),
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) +
    scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
    theme(aspect.ratio=1) +
    ggtitle(title) +
    theme_update(plot.title = element_text(size=16, face="bold", hjust=0.5),
                 axis.title.x = element_text(size=12, face= "bold"),
                 axis.text.x = element_text(size=12, face="bold"),
                 axis.title.y.left = element_text(size=12, face= "bold"),
                 axis.text.y = element_text(size=12, face="bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) 
  
  
  return(volc.plt)
  
}
options(ggrepel.max.overlaps = Inf)
sig_cut <- 0.05
logfc_cut <- 1
titles <- paste0("KOvsWT ", Time)
xrange    <- list(c(-8, 22, 4), c(-8,10,2), c(-4,10,2), c(-10,10, 2)) 
yrange    <- list(c(0, 32, 4), c(0,80,10), c(0,80,10), c(0, 104, 10))
topN      <- 10
xlabel    <- "log2FC (KO/WT)"
ylabel    <- "-log"[10] ~ "(adj.p.value)"
volplt_0w <- functionPlotDEVolcano(RES.files1[1], ensEMBL2id, sig_cut, logfc_cut, titles[1],  xrange[[1]], yrange[[1]], topN, xlabel, ylabel)
volplt_2w <- functionPlotDEVolcano(RES.files1[2], ensEMBL2id, sig_cut, logfc_cut, titles[2],  xrange[[2]], yrange[[2]], topN, xlabel, ylabel)
volplt_4w <- functionPlotDEVolcano(RES.files1[3], ensEMBL2id, sig_cut, logfc_cut, titles[3],  xrange[[3]], yrange[[3]], topN, xlabel, ylabel)
volplt_7w <- functionPlotDEVolcano(RES.files1[4], ensEMBL2id, sig_cut, logfc_cut, titles[4],  xrange[[4]], yrange[[4]], topN, xlabel, ylabel)


pdf("./Results_Tables_Figures/March_2023/CAD_jcg23_0001-VolcanoPlot_0w_lfcshrink_March_2023.pdf")
print(volplt_0w)
dev.off()
pdf("./Results_Tables_Figures/March_2023/CAD_jcg23_0001-VolcanoPlot_2w_lfcshrink_March_2023.pdf")
print(volplt_2w)
dev.off()
pdf("./Results_Tables_Figures/March_2023/CAD_jcg23_0001-VolcanoPlot_4w_lfcshrink_March_2023.pdf")
print(volplt_4w)
dev.off()
pdf("./Results_Tables_Figures/March_2023/CAD_jcg23_0001-VolcanoPlot_7w_lfcshrink_March_2023.pdf")
print(volplt_7w)
dev.off()

message("+------------Selected individual genes plot across all samples ------------+") 

genesname  <- c("Ppp1r15a", "Kcnj11", "Kcnj14", "Ddit3", "Il6", "Gdf15", "Atf4", "Trp53", "Pck1",
                "Asns", "Sesn2", "Trib3", "Xbp1", "Atf3", "Psat1", "Psph", "Phgdh", "Bcat2",
                "Mycn", "Slc7a11", "Gss", "Spic", "Slc1a4", "Nrp1")
genes2plot <- lapply(genesname, function(x) ensEMBL2id[which(ensEMBL2id$external_gene_name==x), "ensembl_gene_id"])
genes2plot <- unlist(genes2plot)[-5]
colData(dds_all)$ComCon <- paste0(colData(dds_all)$Condition, "_", colData(dds_all)$Time)
makeGeneCountPlot <- function(DDS, ensEMBL2id, CONDITION, gene2plot, genename) {
  # Plot the log2 normalised read counts for a specified gene
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2) <- c("count", "condition")
  t2$count <- log2(t2$count)
  t2$condition <- factor(t2$condition, levels=c("WT_0w", "WT_2w","WT_4w","WT_7w",
                                                "KO_0w", "KO_2w","KO_4w","KO_7w"))
  plt.cont <-  ggplot(t2, aes(x=condition, y=count, fill=condition)) + 
    geom_boxplot(width = 0.5, color=rep(c("blue", "red"), each=4), alpha=0.5, outlier.shape=NA) + 
    geom_point(aes(fill=condition),position=position_jitterdodge(),size=0.6, alpha=0.5) +
    scale_fill_manual(name="Condition", values = rep(c("blue", "red"), each=4)) +
    scale_color_manual(name="Condition", values = rep(c("blue", "red"), each=4)) +
    theme(text = element_text(size=elementTextSize), legend.position="none") +
    ggtitle(genename) + 
    xlab("") + ylab("log2(Normalised count)") +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=elementTextSize,angle = 45, hjust = 1))
  
  
  t2$samples   <- rownames(t2)
  colnames(t2) <- c("count", "condition", "samples")
  t2           <- t2[order(t2$condition),]
  t2$samples2 <- factor(t2$samples, as.character(t2$samples))
  
  plt.ind <- ggplot(t2, aes(x=samples2, y=count, fill=condition, group=condition)) + 
    geom_bar(stat="identity", alpha=0.5) +
    scale_fill_manual(name="Condition", values = rep(c("blue", "red"), each=4)) +
    xlab("") + ylab("log2(Normalised count)") +
    ggtitle(genename) +
    theme_classic() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=elementTextSize, angle = 45, hjust = 1))
  
  
  print(paste("Created plot for", gene2plot), sep=" ")
  
  return(list(plt.cont, plt.ind))
  }


pdf("./Results_Tables_Figures/March_2023/CAD_jcg23_0001-SelGene_Indivi_Group_Plot_March_2023.pdf", width=14, height = 14)
ind.plt1 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[1], genesname[1])
ind.plt2 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[2], genesname[2])
ind.plt3 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[3], genesname[3])
ind.plt4 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[4], genesname[4])
ind.plt5 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[5], genesname[5])
ind.plt6 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[6], genesname[6])
ind.plt7 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[7], genesname[7])
ind.plt8 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[8], genesname[8])
ind.plt9 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[9], genesname[9])

ind.plt10 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[10], genesname[10])
ind.plt11 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[11], genesname[11])
ind.plt12 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[12], genesname[12])
ind.plt13 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[13], genesname[13])
ind.plt14 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[14], genesname[14])
ind.plt15 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[15], genesname[15])
ind.plt16 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[16], genesname[16])
ind.plt17 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[17], genesname[17])
ind.plt18 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[18], genesname[18])
ind.plt19 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[19], genesname[19])

ind.plt20 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[20], genesname[20])
ind.plt21 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[21], genesname[21])
ind.plt22 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[22], genesname[22])
ind.plt23 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[23], genesname[23])
ind.plt24 <- makeGeneCountPlot(dds_all, ensEMBL2id, "ComCon", genes2plot[24], genesname[24])

plot_grid(ind.plt1[[1]], ind.plt1[[2]], ind.plt2[[1]], ind.plt2[[2]],
          ind.plt3[[1]], ind.plt3[[2]], ind.plt4[[1]], ind.plt4[[2]],
          nrow=4, byrow=F)
plot_grid(ind.plt5[[1]], ind.plt5[[2]], ind.plt6[[1]], ind.plt6[[2]],
          ind.plt7[[1]], ind.plt7[[2]], ind.plt8[[1]], ind.plt8[[2]],
          nrow=4, byrow=F)

plot_grid(ind.plt9[[1]], ind.plt9[[2]], ind.plt10[[1]], ind.plt10[[2]],
          ind.plt11[[1]], ind.plt11[[2]], ind.plt12[[1]], ind.plt12[[2]],
          nrow=4, byrow=F)
plot_grid(ind.plt13[[1]], ind.plt13[[2]], ind.plt14[[1]], ind.plt14[[2]],
          ind.plt15[[1]], ind.plt15[[2]], ind.plt16[[1]], ind.plt16[[2]],
          nrow=4, byrow=F)
plot_grid(ind.plt17[[1]], ind.plt17[[2]], ind.plt18[[1]], ind.plt18[[2]],
          ind.plt19[[1]], ind.plt19[[2]], ind.plt20[[1]], ind.plt20[[2]],
          nrow=4, byrow=F)
plot_grid(ind.plt21[[1]], ind.plt21[[2]], ind.plt22[[1]], ind.plt22[[2]],
          ind.plt23[[1]], ind.plt23[[2]], ind.plt24[[1]], ind.plt24[[2]],
          nrow=4, byrow=F)

dev.off()


message("+------- Gene Ontology analysis for wk4,7 KOvsWT first using clusterProfiler. modified on 27/02/2023------+")

suppressPackageStartupMessages({
  library("clusterProfiler")
  library("DOSE")
  library("GSEABase")
  library("AnnotationHub")
  library("org.Hs.eg.db")
  library("org.Mm.eg.db")
  library("BiocParallel")
  library("gage")
  library("gageData")
  library("enrichplot")
  library("reshape2")
  library("enrichplot")
  library("ggraph")
  library("ggforce")
  library("reactome.db")
  library("ReactomePA")
  library("Biostrings")
})

hub <- AnnotationHub()
query(hub, "Mus_musculus")
register(MulticoreParam(2))


load("/Users/xz289/Documents/Minoru/Project_mt709_0001/Data/References_Data/GRCm39/Ensembl_GRCm39_entrezID_ensembID_exterName.RData")

## use lfcshrink results
GO_Data_conversion <- function(resfile, ensGO, padj, l2fc1, l2fc2, selCol.names, outfile){
  resdat <- read.csv(resfile, header = T)
  resdat <- subset(resdat, !is.na(padj) & !is.na(log2FoldChange))
  resann <- merge(resdat, ensGO, by = "external_gene_name", all.x = T)
  ## select the first l2fc cutoff
  resann.sig1 <- subset(resann, padj <= padj & abs(log2FoldChange) >= l2fc1)
  resdf1 <- resann.sig1[, selCol.names]
  colnames(resdf1)  <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC")
  resdf1 <- resdf1[order(-resdf1$L2FC),]
  ## select the second l2fc cutoff
  resann.sig2 <- subset(resann, padj <= padj & abs(log2FoldChange) >= l2fc2)
  resdf2 <- resann.sig2[, selCol.names]
  colnames(resdf2)  <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC")
  resdf2 <- resdf2[order(-resdf2$L2FC),]
  
  write.xlsx(resdf1, file = outfile, sheetName = "sigDEGs1", append = F, row.names =F)
  write.xlsx(resdf2, file = outfile, sheetName = "sigDEGs2", append = T, row.names =F)
  
  resdf.bk <- resann[, selCol.names]
  colnames(resdf.bk) <- colnames(resdf1) 
  resdf.bk <- resdf.bk[order(-resdf.bk$L2FC),]
  write.xlsx(resdf.bk, file = outfile, sheetName = "BKDEGs", append = T)
  gc()
}
ensGO <- ensEMBL2id.GO
sel.column.names <-c("entrezgene_id", "ensembl_gene_id.x", "external_gene_name", "log2FoldChange")
outfiles <- paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/GO_input_wk", c(4, 7), "_l2fc0.6_l2fc1_KOvsWT_list-March_2023.xlsx")

wk4_goD <- GO_Data_conversion(RES.files1[3], ensGO, 0.05, 0.6, 1, sel.column.names, outfiles[1])
wk7_goD <- GO_Data_conversion(RES.files1[4], ensGO, 0.05, 0.6, 1, sel.column.names, outfiles[2])

## approach 1, up and down regulated genes all together. Only do BP and KEGG, plus gsea to quantify.
## gene 
## KEGG pathways are not stable, rerun later, 09/03/2032
GO_enrich_fn <- function(input, inputInd, bkInd, pval=0.05, ontology, gsemin=10, gsemax=500, outfile){
  ## input is the data, inputInd is the sheetIndex for input data, bkInd is the index for background.
  ## pval=0.05,  ontology (BP, CC, MF or All)
  ## gsemin= 10, gsemax = 500, outfile(saved R object)
  set.seed(1234)
  resdf       <- read_excel(input, sheet = inputInd)
  resdf.bk    <- read_excel(input, sheet = bkInd)
  geneDat     <- subset(resdf, !is.na(resdf$ENTREZID))
  geneList    <- geneDat$L2FC
  names(geneList) <- geneDat$ENTREZID
  
  geneDat.bk     <- subset(resdf.bk, !is.na(resdf.bk$ENTREZID))
  geneList.bk    <- geneDat.bk$L2FC
  names(geneList.bk) <- geneDat.bk$ENTREZID
  
  enrichGO_Res <- enrichGO(gene = as.character(na.omit(resdf$ENTREZID)), 
                          universe = as.character(na.omit(resdf.bk$ENTREZID)),
                          keyType = "ENTREZID", 
                          OrgDb = org.Mm.eg.db,
                          ont = ontology, 
                          pAdjustMethod = "BH", 
                          pvalueCutoff  = pval,
                          readable      = TRUE)
  
  enrichGSE_Res  <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db,
                  ont = "BP", minGSSize = 10, maxGSSize = 500, eps=0,
                  pvalueCutoff = 0.05, verbose = FALSE, nPermSimple = 10000)
  if(dim(as.data.frame(enrichGSE_Res))[1]!=0){
    enrichGSE_Res  <- setReadable(enrichGSE_Res, OrgDb = org.Mm.eg.db, keyType="ENTREZID")}
  
  enrichKEGG_Res <-  enrichKEGG(gene = names(geneList), universe=names(geneList.bk), organism = 'mmu', pvalueCutoff = 0.05)
  if(dim(as.data.frame(enrichKEGG_Res))[1]!=0){
  enrichKEGG_Res <- setReadable(enrichKEGG_Res, OrgDb = org.Mm.eg.db, keyType="ENTREZID")}
  
  gseKEGG_Res <- gseKEGG(geneList, organism = "mmu", keyType = "kegg",exponent = 1,
                          minGSSize = 10, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05)
  if(dim(as.data.frame(gseKEGG_Res))[1]!=0){
    gseKEGG_Res <- setReadable(gseKEGG_Res, OrgDb = org.Mm.eg.db, keyType="ENTREZID")}
  
  save(enrichGO_Res, enrichGSE_Res, enrichKEGG_Res, gseKEGG_Res, file = outfile)
       #enrichKEGG_Res, 
  gc()
}
gooutfiles1 <- paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/GO_Analysis_wk", c(4, 7), "_l2fc0.6_KOvsWT_BP_KEGG-March_2023.RData")

wk4_GOAnalysis <- GO_enrich_fn(outfiles[1], 1, 3, pval=0.05,  "All",  gsemin=10, gsemax=500, gooutfiles1[1])
wk7_GOAnalysis <- GO_enrich_fn(outfiles[2], 1, 3, pval=0.05,  "All",  gsemin=10, gsemax=500, gooutfiles1[2])

load(gooutfiles1[2])
write.xlsx(as.data.frame(enrichKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk7_Wk4_KOvsWT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=F, sheetName = "enrichKEGG_wk7")
write.xlsx(as.data.frame(gseKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk7_Wk4_KOvsWT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "gseKEGG_wk7")
write.xlsx(as.data.frame(enrichGO_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk7_Wk4_KOvsWT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=F, sheetName = "wk7_enrichGO")
write.xlsx(as.data.frame(enrichGSE_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk7_Wk4_KOvsWT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "wk7_enrichGSE_BP")

load(gooutfiles1[1])
write.xlsx(as.data.frame(enrichKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk7_Wk4_KOvsWT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichKEGG_wk4")
write.xlsx(as.data.frame(enrichGO_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk7_Wk4_KOvsWT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "wk4_enrichGO")

message("+-----semantic similarity among GO terms  rrvgo package (1.8.0), esp BP, enrichGO and enrichGSE-----+")
message("+Reduce and visualize lists of Gene Ontology terms by identifying redudance based on semantic similarity.+")

library("pathview")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("org.Mm.eg.db")
library("rrvgo")
library("GOSemSim")


GOSeSimdata.BP   <- godata(OrgDb = "org.Mm.eg.db", keytype = "ENTREZID", ont="BP", computeIC = TRUE)

GOSemSim_Data_function <- function(goterm, gofun, model, goName, threshcut, Project, max.overlap){
  ego.all       <- goterm
  ego.dat       <- as.data.frame(ego.all)
  ego.dat       <- subset(ego.dat, ONTOLOGY==gofun) 
  ego.dat.dim   <- dim(ego.dat)[1]
  
  simMatrix     <- calculateSimMatrix(ego.dat$ID, orgdb="org.Mm.eg.db", ont=gofun, semdata=GOSeSimdata.BP, method="Rel")
  scores        <- setNames(-log10(ego.dat$qvalue), ego.dat$ID)
  reducedTerms  <- reduceSimMatrix(simMatrix, scores, threshold=threshcut, orgdb="org.Mm.eg.db")
  maxclust      <- max(reducedTerms$cluster)
  
  message("+---------               Supplementary data for semantic similarity                   --------------+")
  
  write.csv(reducedTerms, file = paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/",Project, "-", model, "_GOSeSim_", goName, "_", gofun, "_N", ego.dat.dim, "_thresh_", threshcut, "_cluster", maxclust, "_summaryTable_March_2023.csv"), row.names=F)
  save(reducedTerms, simMatrix, file = paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/",Project, "-", model, "_GOSeSim_", goName, "_", gofun, "_N", ego.dat.dim, "_thresh_", threshcut, "_cluster", maxclust, "_summaryTable_March_2023.RData"))
  
  options(ggrepel.max.overlaps = max.overlap)
  
  message("+---------               Supplementary Fig for semantic similarit                      --------------+")
  
  pdf(paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/",Project, "-", model, "_GOSeSim_", goName, "_", gofun, "_N", ego.dat.dim, "_thresh_", threshcut, "_cluster", maxclust, "_scatterPlot_March_2023.pdf"))
  scatterPlot(simMatrix, reducedTerms, size="score", labelSize = 2)
  dev.off()
}
GOSemSim_Data_function_GSEBP <- function(goterm, gofun, model, goName, threshcut, Project, max.overlap){
  ego.all       <- goterm
  ego.dat       <- as.data.frame(ego.all)
  ego.dat.dim   <- dim(ego.dat)[1]
  
  simMatrix     <- calculateSimMatrix(ego.dat$ID, orgdb="org.Mm.eg.db", ont=gofun, semdata=GOSeSimdata.BP, method="Rel")
  scores        <- setNames(-log10(ego.dat$qvalue), ego.dat$ID)
  reducedTerms  <- reduceSimMatrix(simMatrix, scores, threshold=threshcut, orgdb="org.Mm.eg.db")
  maxclust      <- max(reducedTerms$cluster)
  
  message("+---------               Supplementary data for semantic similarity                   --------------+")
  
  write.csv(reducedTerms, file = paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/",Project, "-", model, "_GOSeSim_", goName, "_", gofun, "_N", ego.dat.dim, "_thresh_", threshcut, "_cluster", maxclust, "_summaryTable_March_2023.csv"), row.names=F)
  save(reducedTerms, simMatrix, file = paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/",Project, "-", model, "_GOSeSim_", goName, "_", gofun, "_N", ego.dat.dim, "_thresh_", threshcut, "_cluster", maxclust, "_summaryTable_March_2023.RData"))
  
  options(ggrepel.max.overlaps = max.overlap)
  
  message("+---------               Supplementary Fig for semantic similarit                      --------------+")
  
  pdf(paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/",Project, "-", model, "_GOSeSim_", goName, "_", gofun, "_N", ego.dat.dim, "_thresh_", threshcut, "_cluster", maxclust, "_scatterPlot_March_2023.pdf"))
  scatterPlot(simMatrix, reducedTerms, size="score", labelSize = 2)
  dev.off()
}

load(gooutfiles1[1])
allBP.wk4  <- GOSemSim_Data_function(enrichGO_Res,"BP", "wk4", "allBP", 0.9, Project, 160)

load(gooutfiles1[2])
allBP.wk7  <- GOSemSim_Data_function(enrichGO_Res,"BP", "wk7", "allBP", 0.9, Project, 160)
gseBP.wk7  <- GOSemSim_Data_function_GSEBP(enrichGSE_Res, "BP", "wk7", "gseBP", 0.9, Project, 160)

message("+------------reactome only with up/down together-------------------------+")

## prepare data for comparison later stage wk4 and wk7 KOvsWT, up/dw data frame
mydf1 <- read_excel(outfiles[1], sheet = 1); mydf2 <- read_excel(outfiles[2], sheet = 1)
mydf1 <- unique(subset(mydf1, !is.na(ENTREZID)))
mydf2 <- unique(subset(mydf2, !is.na(ENTREZID)))  
mydf1$group <- ifelse(mydf1$L2FC > 0, "upregulated", "downregulated")
mydf1$othergroup <- "wk4"
mydf2$group <- ifelse(mydf2$L2FC > 0, "upregulated", "downregulated")
mydf2$othergroup <- "wk7"

cmydf <- rbind(mydf1, mydf2)
table(cmydf$group, cmydf$othergroup)
##wk4  wk7
##downregulated   20 2625
##upregulated    136 3154
## the final analysis may have slightly number difference due to the missingness of geneID,Name or others.

write.csv(cmydf, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Input_data_compareCluster_wk4swk7_l2fc06_list_March_2023.csv", row.names=F)

## Data conversion for reactome from mouse to human orthology
## get the cmydf geneSymbol to gprofiler2 orthology and convert to human gene name

ensembl_Ortho_HM <- read.delim("./Results_Tables_Figures/March_2023/Gene_Ontology/Human_Mouse_orthology_Ensemble_Biomart_export_March_2023.txt", header=T)
## check the genes in ortho.hm which did not have orhtology gene in ensembl
ortho.hm <- read.csv("./Results_Tables_Figures/March_2023/Gene_Ontology/KOvsWT_wk4_wk7_gProfiler_mmusculus_hsapiens_01-03-2023_16-26-28.csv", header=T)
test1 <- ortho.hm$ortholog_name[-which(ortho.hm$ortholog_name=="N/A")]
test <- test1[which(test1%in%ensembl_Ortho_HM$Human.gene.name==F)]
write.table(as.data.frame(unique(test)), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/testOrtho_list.txt", quote=F, row.names=F, col.names=F)
ensembl_Ortho_HM_1 <- ensembl_Ortho_HM[,-2]
ensembl_Ortho_HM_2 <- read.csv("./Results_Tables_Figures/March_2023/Gene_Ontology/gProfiler_hsapiens_mmusculus_13-03-2023_10-23-01.csv")
ensembl_Ortho_HM_2 <- ensembl_Ortho_HM_2[,c(6,5,3,2)]
colnames(ensembl_Ortho_HM_2) <- colnames(ensembl_Ortho_HM_1)

new_HM_Ortho <- unique(rbind(ensembl_Ortho_HM_1, ensembl_Ortho_HM_2))

write.csv(new_HM_Ortho, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Ensembl_Mouse_Human_FullOrthology_list_March_2023.csv", quote=F, row.names=F)

message("+----Check the overlap for wk4 or wk7 all genes orthology with human name--------+")
message("+----And then convert to entrezeID for reactome analysis                 --------+")

## wk4 data sorting
wk4.df <- read_excel(outfiles[1], sheet = 1); 
wk4.df.bk <- read_excel(outfiles[1], sheet = 3);
colnames(new_HM_Ortho)[1] <- c("ENSEMBL")
wk4.orthos    <- merge(wk4.df, new_HM_Ortho, by = "ENSEMBL"); 
wk4.bk.orthos <- merge(wk4.df.bk, new_HM_Ortho, by = "ENSEMBL");
wk4.orthos    <- wk4.orthos[order(-wk4.orthos$L2FC),]
wk4.bk.orthos <- wk4.bk.orthos[order(-wk4.bk.orthos$L2FC),]
reactome.wk4 <- bitr(unique(wk4.orthos$Human.gene.name), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
reactome.wk4.universe <- bitr(unique(wk4.bk.orthos$Human.gene.name), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

## wk7 data sorting
wk7.df <- read_excel(outfiles[2], sheet = 1)
wk7.df.bk <- read_excel(outfiles[2], sheet = 3);
wk7.orthos    <- merge(wk7.df, new_HM_Ortho, by = "ENSEMBL"); 
wk7.bk.orthos <- merge(wk7.df.bk, new_HM_Ortho, by = "ENSEMBL");
wk7.orthos    <- wk7.orthos[order(-wk7.orthos$L2FC),]
wk7.bk.orthos <- wk7.bk.orthos[order(-wk7.bk.orthos$L2FC),]
reactome.wk7 <- bitr(unique(wk7.orthos$Human.gene.name), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
reactome.wk7.universe <- bitr(unique(wk7.bk.orthos$Human.gene.name), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

message("+---enrichPathway, reactome for up/dw wk4/7 analysis ----------------+")
## BiocManager::install("ReactomePA")
## need to convert mouse gene to human gene.
library("ReactomePA")

wk4_reactome <- enrichPathway(reactome.wk4$ENTREZID,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.2,
                              universe=reactome.wk4.universe$ENTREZID,
                              minGSSize = 10,
                              maxGSSize = 500,
                              readable = TRUE)

wk7_reactome <- enrichPathway(reactome.wk7$ENTREZID,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.2,
                              universe=reactome.wk7.universe$ENTREZID,
                              minGSSize = 10,
                              maxGSSize = 500,
                              readable = TRUE)

message("+---GSE pathways of reactome analysis for wk4 and wk7----------+")

colnames(reactome.wk4)[1] <- "Human.gene.name"
wk4_geneList <- merge(wk4.orthos, reactome.wk4, by = "Human.gene.name")
wk4_geneList <- wk4_geneList[order(-wk4_geneList$L2FC),]
wk4_rgeneList <- wk4_geneList$L2FC
names(wk4_rgeneList) <- wk4_geneList$ENTREZID.y
wk4_rgeneList <- wk4_rgeneList[!duplicated(names(wk4_rgeneList))]
wk4_gseReactome <- gsePathway(wk4_rgeneList, 
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE) ## no enrichment even I changed to 0.2

colnames(reactome.wk7)[1] <- "Human.gene.name"
wk7_geneList <- merge(wk7.orthos, reactome.wk7, by = "Human.gene.name")
wk7_geneList <- wk7_geneList[order(-wk7_geneList$L2FC),]
wk7_rgeneList <- wk7_geneList$L2FC
names(wk7_rgeneList) <- wk7_geneList$ENTREZID.y
wk7_rgeneList <- wk7_rgeneList[!duplicated(names(wk7_rgeneList))]
wk7_gseReactome <- gsePathway(wk7_rgeneList, 
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH", 
                              minGSSize = 10, 
                              maxGSSize = 500, 
                              verbose = FALSE)
wk7_gseReactome <- setReadable(wk7_gseReactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.xlsx(as.data.frame(wk4_reactome), file= "./Results_Tables_Figures/March_2023/Gene_Ontology/Reactome_updw_enrich_gse_wk4_wk7_summaryTable.xlsx",
           append=F, row.names=F, sheetName = "wk4_enrichPathway")
write.xlsx(as.data.frame(wk7_reactome), file= "./Results_Tables_Figures/March_2023/Gene_Ontology/Reactome_updw_enrich_gse_wk4_wk7_summaryTable.xlsx",
           append=T, row.names=F, sheetName = "wk7_enrichPathway")
write.xlsx(as.data.frame(wk7_gseReactome), file= "./Results_Tables_Figures/March_2023/Gene_Ontology/Reactome_updw_enrich_gse_wk4_wk7_summaryTable.xlsx",
           append=T, row.names=F, sheetName = "wk7_gsePathway_new")
save(wk4_reactome, wk7_reactome, wk7_gseReactome, wk7_rgeneList, wk4_rgeneList,
     reactome.wk7, reactome.wk4, reactome.wk7.universe, reactome.wk4.universe,
     file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Reactome_updw_enrich_gse_wk4_wk7_summaryTable.RData")
save(wk4_geneList, wk7_geneList,
     file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Reactome_inputortho_KO_wk4_wk7_data.RData")


message("+----- Compare clusters gene ontology with BP and Reactome for two approaches------+")
message("+----Only compare wk4 and wk7, without split up and dw regulated genes ------------+")

formula_GOAll <- compareCluster(ENTREZID~othergroup, data=cmydf, fun="enrichGO", OrgDb = org.Mm.eg.db, ont = "BP")
formula_GOAll <- setReadable(formula_GOAll, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

## mouse ID convert to human ID, gProfiler2. 
my.symbols <- na.omit(ortho.hm$ortholog_name)
reactID <- bitr(my.symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
## merge back to mouse data and groups
colnames(ortho.hm)[c(3,5,6)] <- c("ENSEMBL", "SYMBOL_H", "ENSEMBL_H")
colnames(reactID) <- c("SYMBOL_H", "ENTREZID_H")
merreact1 <- merge(ortho.hm[, c(2,3,5,6)], reactID, by = "SYMBOL_H")
merreact2 <- merge(cmydf, merreact1, by = "ENSEMBL")
merreact2$mgroup <- paste0(merreact2$othergroup,"_", merreact2$group)
merreact31 <- merreact2 %>% arrange(othergroup, -L2FC)
merreact32 <- merreact2 %>% arrange(mgroup, -L2FC)
formula_ReactAll <- compareCluster(ENTREZID_H~othergroup, data=merreact31, fun="enrichPathway", readable=TRUE)
## wk4 has immune and cellular response to stimuli
## wk7has immune, metabolism main, and extracellualr matrix organisation(EMO), a few immune, protein localisation
## metabolism and protein, disease, hemostasis, et al.

write.xlsx(as.data.frame(formula_GOAll), 
           file = "./Results_Tables_Figures/March_2023/Gene_Ontology/CompareCluster_enrichGO_pathway_wk4_wk7_KOvsWT_all_list_01_03_2023.xlsx", row.names=F,
           append=F, sheetName = "EnrichGO_compare_all")
write.xlsx(as.data.frame(formula_ReactAll), 
           file = "./Results_Tables_Figures/March_2023/Gene_Ontology/CompareCluster_Reactome_pathway_wk4_wk7_KOvsWT_all_list_01_03_2023.xlsx", row.names=F,
           append=F, sheetName = "Reactome_compare_all")


message("+-----BP, reactome compare for wk4 and wk7 up/dw separately---------------+")

formula_GO <- compareCluster(ENTREZID~group+othergroup, data=cmydf, fun="enrichGO", OrgDb = org.Mm.eg.db, ont = "BP")
formula_GO <- setReadable(formula_GO, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

write.xlsx(as.data.frame(formula_GO), 
           file = "./Results_Tables_Figures/March_2023/Gene_Ontology/CompareCluster_enrichGO_pathway_wk4_wk7_KOvsWT_all_list_01_03_2023.xlsx", row.names=F,
           append=T, sheetName = "EnrichGO_compare_up_down1")
message("+---- Reactom compare pathway analysis for wk4/7 up/dw regulated genes, l2fc=0.6 -----------------+")

formula_React <- compareCluster(ENTREZID_H~group+othergroup, data=merreact32, fun="enrichPathway", readable=TRUE)
write.xlsx(as.data.frame(formula_React), 
           file = "./Results_Tables_Figures/March_2023/Gene_Ontology/CompareCluster_Reactome_pathway_wk4_wk7_KOvsWT_all_list_01_03_2023.xlsx",
           row.names=F, append=T, sheetName = "Reactome_compare_up_down")

pdf("./Results_Tables_Figures/Reactome_dotplot_wk4_wk7_up_dw_pathways.pdf", height =100, width =12)
dotplot(formula_React, showCategory=150)
dev.off()

message("+-----KEGG compare for wk4 and wk7 up/dw separately---------------+")

formula_KEGGAll <- compareCluster(ENTREZID~othergroup, data=cmydf, fun="enrichKEGG", organism = "mmu")
formula_KEGGAll <- setReadable(formula_KEGGAll, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

formula_KEGGsplit <- compareCluster(ENTREZID~group+othergroup, data=cmydf, fun="enrichKEGG", organism = "mmu")
formula_KEGGsplit <- setReadable(formula_KEGGsplit, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

write.xlsx(as.data.frame(formula_KEGGAll), 
           file = "./Results_Tables_Figures/March_2023/Gene_Ontology/CompareCluster_enrichGO_pathway_wk4_wk7_KOvsWT_all_list_01_03_2023.xlsx", row.names=F,
           append=T, sheetName = "EnrichKEGG_compare")
write.xlsx(as.data.frame(formula_KEGGsplit), 
           file = "./Results_Tables_Figures/March_2023/Gene_Ontology/CompareCluster_enrichGO_pathway_wk4_wk7_KOvsWT_all_list_01_03_2023.xlsx", row.names=F,
           append=T, sheetName = "EnrichKEGG_compare_up_down")
save(formula_KEGGAll, formula_KEGGsplit, formula_React, formula_GO, formula_ReactAll,formula_GOAll, merreact31,merreact32,
     file = "./Results_Tables_Figures/March_2023/Gene_Ontology/CompareCluster_wk47_KOvsWT_all_splitupdw_April_2023.RData")

message("+-------------------------------------------------------------------------------------------------+")
message("+---- check across WT time points analysis                                       -----------------+")
message("+-------------------------------------------------------------------------------------------------+")
customPCA_pairT     <- function(sampleTBL, ensEMBL2id, RLD, TOPNUM, GRPName) {
  RLD            <- as.data.frame(RLD) 
  RLD$ensembl_gene_id <- rownames(RLD)
  RLD.mer<- merge(RLD, ensEMBL2id, by = "ensembl_gene_id")
  RLD.mer<- RLD.mer[-which(duplicated(RLD.mer$external_gene_name)==T),]
  rownames(RLD.mer) <- RLD.mer$external_gene_name
  RLD.mer<- RLD.mer[,grep(GRPName, colnames(RLD.mer))]
  rv     <- rowVars(as.matrix(RLD.mer))
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD.mer[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sample, pca$x, 
                          condition=sampleTBL$Time)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3, alpha=0.75) + 
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
    theme(text = element_text(size=elementTextSize)) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3, alpha=0.75 ) + 
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(colour=condition, fill=condition, group = condition)) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize))+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  loadings                 <- as.data.frame(pca$rotation)
  pca.1         <- loadings[order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:25),]
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(rownames(pca.1.25),levels=unique(rownames(pca.1.25))), y=PC1)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  pca.2         <- loadings[order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:25),]
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(rownames(pca.2.25),levels=unique(rownames(pca.2.25))), y=PC2)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme()+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  return(list(plt.pca,plt.pca.nl,pca.1.25.plot, pca.2.25.plot))
  
}
Time1 <- c("0w", "0w", "0w", "2w", "2w", "4w")
Time2 <- c("2w", "4w", "7w", "4w", "7w", "7w")

RESsigW.files <- paste0("./Results_Tables_Figures/March_2023/", Project, "-ResSig_WT_", Time2, "vs", Time1, "_summaryTable_March_2023.csv")
RESW.files <- paste0("./Results_Tables_Figures/March_2023/",Project, "-Resall_WT_",Time2, "vs", Time1, "_summaryTable_March_2023.csv")
RDW.files  <- paste0("./Results_Tables_Figures/March_2023/",Project, "-RData_WT_",Time2, "vs", Time1,"_March_2023.RData")

RESsigW1.files <- paste0("./Results_Tables_Figures/March_2023/",Project, "-ResSig_WT_lfcshrink_", Time2, "vs", Time1, "_summaryTable_March_2023.csv")
RESW1.files <- paste0("./Results_Tables_Figures/March_2023/",Project, "-Resall_WT_lfcshrink_",Time2, "vs", Time1, "_summaryTable_March_2023.csv")
RDW1.files  <- paste0("./Results_Tables_Figures/March_2023/",Project, "-RData_WT_lfcshrink_",Time2, "vs", Time1,"_March_2023.RData")

count_df_allWT <- count_df_all[,grep("WT", colnames(count_df_all))]
targets_allWT  <- targets_all[targets_all$Condition=="WT",]

for(i in 1:length(Time1)){
  grp      <- paste0("T", Time1[i], "|T", Time2[i])
  count_df <- count_df_allWT[,grep(grp, colnames(count_df_allWT))]
  targets_sub  <- targets_allWT[grep(grp,targets_allWT$sample),]
  dds_sub   <- DESeqDataSetFromMatrix(countData=count_df, colData=targets_sub, design=~Time)
  dds_sub$Time <- relevel(dds_sub$Time, ref = Time1[i])
  dds_sub   <- DESeq(dds_sub, parallel=TRUE)
  dds_sub   <- estimateSizeFactors(dds_sub)
  vsd_sub   <- vst(dds_sub,     blind=F)
  colData(vsd_sub)
  resultsNames(dds_sub)
  normalized_counts <- as.data.frame(counts(dds_sub, normalized=TRUE))
  normalized_counts$ensembl_gene_id <- rownames(normalized_counts)
  normcounts <- merge(normalized_counts, ensEMBL2id, by ="ensembl_gene_id")
  
  pdf(paste0("./Results_Tables_Figures/March_2023/CAD_jcg23_0001-DESeq2_pairwisePCA_WT_", Time2[i], "vs", Time1[i], "_March_2023.pdf"), width=8, height= 6)
  pca <- customPCA_pairT(targets_sub, ensEMBL2id, assay(vsd_sub), TOPNUM, "WT")
  plot_grid(pca[[1]], pca[[2]], pca[[3]], pca[[4]], nrow=2, byrow=T)
  dev.off()
  res0       <- results(dds_sub, contrast = c("Time", Time2[i], Time1[i]))
  res0.sig   <- as.data.frame(subset(res0, padj <= significance & abs(log2FoldChange) >= l2fc))
  print(nrow(res0.sig)); print(table(sign(res0.sig$log2FoldChange)))
  test       <- as.data.frame(subset(res0, padj <= significance & abs(log2FoldChange) >= 0.6))
  print(nrow(test)); print(table(sign(test$log2FoldChange)))
  res0.sig$ensembl_gene_id <- rownames(res0.sig)
  res0.sigmer <- merge(res0.sig, ensEMBL2id, by ="ensembl_gene_id")
  write.csv(res0.sigmer, file = RESsigW.files[i])
  res0.merEC <- merECRes_function_pair(res0, ensEMBL2id, dds_sub, RESW.files[i])
  
  res       <- lfcShrink(dds_sub, coef=2, type = "apeglm")
  res.sig   <- as.data.frame(subset(res, padj <= significance & abs(log2FoldChange) >= l2fc))
  print(nrow(res.sig)); print(table(sign(res.sig$log2FoldChange)))
  test1     <- as.data.frame(subset(res, padj <= significance & abs(log2FoldChange) >= 0.6))
  print(nrow(test1)); print(table(sign(test1$log2FoldChange)))
  res.sig$ensembl_gene_id <- rownames(res.sig)
  res.sigmer <- merge(res.sig, ensEMBL2id, by ="ensembl_gene_id")
  write.csv(res.sigmer, file = RESsigW1.files[i])
  res.merEC <- merECRes_function_pair(res, ensEMBL2id, dds_sub, RESW1.files[i])
  save(dds_sub, vsd_sub, normcounts, file=RDW1.files[i])
}

## generate volcano plot for WT 2/4/7 vs 0w
options(ggrepel.max.overlaps = Inf)
sig_cut <- 0.05
logfc_cut <- 1
titles <- c("WT 2wvs0w", "WT 4wvs0w", "WT 7wvs0w")
xrange    <- list(c(-8, 8, 4), c(-8,8,2), c(-6,6,2)) 
yrange    <- list(c(0, 28, 4), c(0,50,10), c(0,80,10))
topN      <- 10
xlabel    <- c("log2FC (2w/0w)","log2FC (4w/0w)", "log2FC (7w/0w)")
ylabel    <- "-log"[10] ~ "(adj.p.value)"
volplt_2wvs0w <- functionPlotDEVolcano(RESW1.files[1], ensEMBL2id, sig_cut, logfc_cut, titles[1],  xrange[[1]], yrange[[1]], topN, xlabel[1], ylabel)
volplt_4wvs0w <- functionPlotDEVolcano(RESW1.files[2], ensEMBL2id, sig_cut, logfc_cut, titles[2],  xrange[[2]], yrange[[2]], topN, xlabel[2], ylabel)
volplt_7wvs0w <- functionPlotDEVolcano(RESW1.files[3], ensEMBL2id, sig_cut, logfc_cut, titles[3],  xrange[[3]], yrange[[3]], topN, xlabel[3], ylabel)
pdf("./Results_Tables_Figures/March_2023/WT_volcanoPlot_w247vsw0_Feb_2022.pdf")
plot(volplt_2wvs0w);
plot(volplt_4wvs0w);
plot(volplt_7wvs0w);
dev.off()


message("+---- DESeq2 analysis based on wild type across the time-------------------------+")

resultsNames(dds_WT)
resWT.2wvs0w      <- lfcShrink(dds_WT, coef="Time_2w_vs_0w", type="apeglm")
resWT.2wvs0w.sig  <- as.data.frame(subset(resWT.2wvs0w, padj <= significance & abs(log2FoldChange) >= l2fc))
print(nrow(resWT.2wvs0w.sig)); print(table(sign(resWT.2wvs0w.sig$log2FoldChange)))
resWT.2wvs0w.sig$ensembl_gene_id <- rownames(resWT.2wvs0w.sig)
resWT.2wvs0w.sigmer <- merge(resWT.2wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
print(nrow(resWT.2wvs0w.sig)); print(table(sign(resWT.2wvs0w.sig$log2FoldChange)))
resWT.2wvs0w.sig$ensembl_gene_id <- rownames(resWT.2wvs0w.sig)
resWT.2wvs0w.sigmer <- merge(resWT.2wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
write.xlsx(resWT.2wvs0w.sigmer, file = "./Results_Tables_Figures/March_2023/WT_timemodel_wk247vs0_sigDEGs_l2fc1_Mar_2023.xlsx",
           row.names=F, append = F, sheetName= "WT_2w_vs_0w_sig")
RESWT.files <- paste0("./Results_Tables_Figures/March_2023/WT_timemodel_wk", c(2,4, 7), "vs0_allDEGs_Mar_2023.csv")
resWT.2wvs0w.merEC <- merECRes_function_pair(resWT.2wvs0w, ensEMBL2id, dds_WT, RESWT.files[1])

##
resWT.4wvs0w      <- lfcShrink(dds_WT, coef="Time_4w_vs_0w", type="apeglm")
resWT.4wvs0w.sig  <- as.data.frame(subset(resWT.4wvs0w, padj <= significance & abs(log2FoldChange) >= l2fc))
print(nrow(resWT.4wvs0w.sig)); print(table(sign(resWT.4wvs0w.sig$log2FoldChange)))
resWT.4wvs0w.sig$ensembl_gene_id <- rownames(resWT.4wvs0w.sig)
resWT.4wvs0w.sigmer <- merge(resWT.4wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
print(nrow(resWT.4wvs0w.sig)); print(table(sign(resWT.4wvs0w.sig$log2FoldChange)))
resWT.4wvs0w.sig$ensembl_gene_id <- rownames(resWT.4wvs0w.sig)
resWT.4wvs0w.sigmer <- merge(resWT.4wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
write.xlsx(resWT.4wvs0w.sigmer, file = "./Results_Tables_Figures/March_2023/WT_timemodel_wk247vs0_sigDEGs_l2fc1_Mar_2023.xlsx",
           row.names=F, append = T, sheetName= "WT_4w_vs_0w_sig")
resWT.4wvs0w.merEC <- merECRes_function_pair(resWT.4wvs0w, ensEMBL2id, dds_WT, RESWT.files[2])

##
resWT.7wvs0w      <- lfcShrink(dds_WT, coef="Time_7w_vs_0w", type="apeglm")
resWT.7wvs0w.sig  <- as.data.frame(subset(resWT.7wvs0w, padj <= significance & abs(log2FoldChange) >= l2fc))
print(nrow(resWT.7wvs0w.sig)); print(table(sign(resWT.7wvs0w.sig$log2FoldChange)))
resWT.7wvs0w.sig$ensembl_gene_id <- rownames(resWT.7wvs0w.sig)
resWT.7wvs0w.sigmer <- merge(resWT.7wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
print(nrow(resWT.7wvs0w.sig)); print(table(sign(resWT.7wvs0w.sig$log2FoldChange)))
resWT.7wvs0w.sig$ensembl_gene_id <- rownames(resWT.7wvs0w.sig)
resWT.7wvs0w.sigmer <- merge(resWT.7wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
write.xlsx(resWT.7wvs0w.sigmer, file = "./Results_Tables_Figures/March_2023/WT_timemodel_wk247vs0_sigDEGs_l2fc1_Mar_2023.xlsx",
           row.names=F, append = T, sheetName= "WT_7w_vs_0w_sig")
resWT.7wvs0w.merEC <- merECRes_function_pair(resWT.7wvs0w, ensEMBL2id, dds_WT, RESWT.files[3])


message("+---- GO compare for wk2,4,7 vs wk0 for BP and KEGG and Reactom------------------+")

outfilesWT <- paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/GO_input_wk", c(2,4, 7), "vswk0_l2fc0.6_l2fc1_WT_list-March_2023.xlsx")
wk2_goDWT <- GO_Data_conversion(RESWT.files[1], ensGO, 0.05, 0.6, 1, sel.column.names, outfilesWT[1])
wk4_goDWT <- GO_Data_conversion(RESWT.files[2], ensGO, 0.05, 0.6, 1, sel.column.names, outfilesWT[2])
wk7_goDWT <- GO_Data_conversion(RESWT.files[3], ensGO, 0.05, 0.6, 1, sel.column.names, outfilesWT[3])

GO_enrich_fn <- function(input, inputInd, bkInd, pval=0.05, ontology, gsemin=10, gsemax=500, outfile){
  ## input is the data, inputInd is the sheetIndex for input data, bkInd is the index for background.
  ## pval=0.05,  ontology (BP, CC, MF or All)
  ## gsemin= 10, gsemax = 500, outfile(saved R object)
  set.seed(1234)
  resdf       <- read_excel(input, sheet = inputInd)
  resdf.bk    <- read_excel(input, sheet = bkInd)
  geneDat     <- subset(resdf, !is.na(resdf$ENTREZID))
  geneList    <- geneDat$L2FC
  names(geneList) <- geneDat$ENTREZID
  
  geneDat.bk     <- subset(resdf.bk, !is.na(resdf.bk$ENTREZID))
  geneList.bk    <- geneDat.bk$L2FC
  names(geneList.bk) <- geneDat.bk$ENTREZID
  
  enrichGO_Res <- enrichGO(gene = as.character(na.omit(resdf$ENTREZID)), 
                           universe = as.character(na.omit(resdf.bk$ENTREZID)),
                           keyType = "ENTREZID", 
                           OrgDb = org.Mm.eg.db,
                           ont = ontology, 
                           pAdjustMethod = "BH", 
                           pvalueCutoff  = pval,
                           readable      = TRUE)
  
  enrichGSE_Res  <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db,
                          ont = "BP", minGSSize = 10, maxGSSize = 500, eps=0,
                          pvalueCutoff = 0.05, verbose = FALSE, nPermSimple = 10000)
  if(dim(as.data.frame(enrichGSE_Res))[1]!=0){
    enrichGSE_Res  <- setReadable(enrichGSE_Res, OrgDb = org.Mm.eg.db, keyType="ENTREZID")}
  
  enrichKEGG_Res <-  enrichKEGG(gene = names(geneList), universe=names(geneList.bk), organism = 'mmu', pvalueCutoff = 0.05)
  if(dim(as.data.frame(enrichKEGG_Res))[1]!=0){
    enrichKEGG_Res <- setReadable(enrichKEGG_Res, OrgDb = org.Mm.eg.db, keyType="ENTREZID")}
  
  gseKEGG_Res <- gseKEGG(geneList, organism = "mmu", keyType = "kegg",exponent = 1,
                         minGSSize = 10, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05)
  if(dim(as.data.frame(gseKEGG_Res))[1]!=0){
    gseKEGG_Res <- setReadable(gseKEGG_Res, OrgDb = org.Mm.eg.db, keyType="ENTREZID")}
  
  save(enrichGO_Res, enrichGSE_Res, enrichKEGG_Res, gseKEGG_Res, file = outfile)
  
  gc()
}
gooutfilesWT <- paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/GO_Analysis_wk", c(2,4, 7), "vswk0_l2fc0.6_WT_BP_KEGG-March_2023.RData")

wk2vs0_GOAnalysis <- GO_enrich_fn(outfilesWT[1], 1, 3, pval=0.05,  "All",  gsemin=10, gsemax=500, gooutfilesWT[1])
wk4vs0_GOAnalysis <- GO_enrich_fn(outfilesWT[2], 1, 3, pval=0.05,  "All",  gsemin=10, gsemax=500, gooutfilesWT[2])
wk7vs0_GOAnalysis <- GO_enrich_fn(outfilesWT[3], 1, 3, pval=0.05,  "All",  gsemin=10, gsemax=500, gooutfilesWT[3])

gooutfilesWT1 <- paste0("./Results_Tables_Figures/March_2023/Gene_Ontology/GO_Analysis_wk", c(2,4, 7), "vswk0_l2fc1_WT_BP_KEGG-March_2023.RData")

wk2vs0_GOAnalysis1 <- GO_enrich_fn(outfilesWT[1], 2, 3, pval=0.05,  "All",  gsemin=10, gsemax=500, gooutfilesWT1[1])
wk4vs0_GOAnalysis1 <- GO_enrich_fn(outfilesWT[2], 2, 3, pval=0.05,  "All",  gsemin=10, gsemax=500, gooutfilesWT1[2])
wk7vs0_GOAnalysis1 <- GO_enrich_fn(outfilesWT[3], 2, 3, pval=0.05,  "All",  gsemin=10, gsemax=500, gooutfilesWT1[3])

load(gooutfilesWT[1])
write.xlsx(as.data.frame(enrichKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=F, sheetName = "enrichKEGG_wk2vswk0_l2fc06")
write.xlsx(as.data.frame(gseKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "gseKEGG_wk2vswk0_l2fc06")

write.xlsx(as.data.frame(enrichGO_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=F, sheetName = "enrichGO_wk2vswk0_l2fc06")
write.xlsx(as.data.frame(enrichGSE_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGSE_BP_wk2vswk0_l2fc06")

load(gooutfilesWT1[1])
write.xlsx(as.data.frame(enrichKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichKEGG_wk2vswk0_l2fc1")
write.xlsx(as.data.frame(gseKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "gseKEGG_wk2vswk0_l2fc1")

write.xlsx(as.data.frame(enrichGO_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGO_wk2vswk0_l2fc1")
#write.xlsx(as.data.frame(enrichGSE_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGSE_BP_wk2vswk0_l2fc1")

## 
load(gooutfilesWT[2])
write.xlsx(as.data.frame(enrichKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichKEGG_wk4vswk0_l2fc06")
write.xlsx(as.data.frame(gseKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "gseKEGG_wk4vswk0_l2fc06")

write.xlsx(as.data.frame(enrichGO_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGO_wk4vswk0_l2fc06")
#write.xlsx(as.data.frame(enrichGSE_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGSE_BP_wk4vswk0_l2fc06")

load(gooutfilesWT1[2])
write.xlsx(as.data.frame(enrichKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichKEGG_wk4vswk0_l2fc1")
#write.xlsx(as.data.frame(gseKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "gseKEGG_wk4vswk0_l2fc1")

write.xlsx(as.data.frame(enrichGO_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGO_wk4vswk0_l2fc1")
write.xlsx(as.data.frame(enrichGSE_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGSE_BP_wk4vswk0_l2fc1")

##
load(gooutfilesWT[3])
write.xlsx(as.data.frame(enrichKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichKEGG_wk7vswk0_l2fc06")
write.xlsx(as.data.frame(gseKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "gseKEGG_wk7vswk0_l2fc06")

write.xlsx(as.data.frame(enrichGO_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGO_wk7vswk0_l2fc06")
write.xlsx(as.data.frame(enrichGSE_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGSE_BP_wk7vswk0_l2fc06")

load(gooutfilesWT1[3])
write.xlsx(as.data.frame(enrichKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichKEGG_wk7vswk0_l2fc1")
write.xlsx(as.data.frame(gseKEGG_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_KEGG_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "gseKEGG_wk7vswk0_l2fc1")

write.xlsx(as.data.frame(enrichGO_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGO_wk7vswk0_l2fc1")
write.xlsx(as.data.frame(enrichGSE_Res), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Wk247vsWk0_WT_GO_enrich_gse_March_2023.xlsx", row.names=F, append=T, sheetName = "enrichGSE_BP_wk7vswk0_l2fc1")

message("+----Reactome for all wk2/4/7 pathway with up/dw together                  -----------------+")
## wk2,4,7vswk0 data sorting
WT.resdf <- list(); WT.resdf.bk <- list();
WT.geneList <- list(); WT.orthos <- list();


for(i in 1:3){
  resdf <- read_excel(outfilesWT[i], sheet = 1); 
  resdf.bk <- read_excel(outfilesWT[i], sheet = 3);
  wt.orthos    <- merge(resdf, new_HM_Ortho, by = "ENSEMBL"); 
  wt.bk.orthos <- merge(resdf.bk, new_HM_Ortho, by = "ENSEMBL");
  wt.orthos    <- wt.orthos[order(-wt.orthos$L2FC),]
  WT.orthos[[i]] <- wt.orthos
  wt.bk.orthos <- wt.bk.orthos[order(-wt.bk.orthos$L2FC),]
  reactome.wt  <- bitr(unique(wt.orthos$Human.gene.name), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  reactome.wt.universe <- bitr(unique(wt.bk.orthos$Human.gene.name), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  WT.resdf[[i]] <- reactome.wt
  WT.resdf.bk[[i]] <- reactome.wt.universe
  
  colnames(reactome.wt)[1] <- "Human.gene.name"
  wt_geneList <- merge(wt.orthos, reactome.wt, by = "Human.gene.name")
  wt_geneList <- wt_geneList[order(-wt_geneList$L2FC),]
  wt_rgeneList <- wt_geneList$L2FC
  names(wt_rgeneList) <- wt_geneList$ENTREZID.y
  wt_rgeneList <- wt_rgeneList[!duplicated(names(wt_rgeneList))]
  WT.geneList[[i]] <- wt_rgeneList
}

save(WT.resdf, WT.resdf.bk,WT.geneList, WT.orthos, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_reactome_input_bkg_geneList_March_2023.RData" )
wt2_reactome <- enrichPathway(WT.resdf[[1]]$ENTREZID,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.2,
                              universe=WT.resdf.bk[[1]]$ENTREZID,
                              minGSSize = 10,
                              maxGSSize = 500,
                              readable = TRUE)

wt4_reactome <- enrichPathway(WT.resdf[[2]]$ENTREZID,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.2,
                              universe=WT.resdf.bk[[2]]$ENTREZID,
                              minGSSize = 10,
                              maxGSSize = 500,
                              readable = TRUE)

wt7_reactome <- enrichPathway(WT.resdf[[3]]$ENTREZID,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.2,
                              universe=WT.resdf.bk[[3]]$ENTREZID,
                              minGSSize = 10,
                              maxGSSize = 500,
                              readable = TRUE)


wt2_gseReactome <- gsePathway(WT.geneList[[1]], 
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH", 
                              verbose = FALSE) 
wt2_gseReactome <- setReadable(wt2_gseReactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

wt4_gseReactome <- gsePathway(WT.geneList[[2]], 
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH", 
                              verbose = FALSE) 
wt4_gseReactome <- setReadable(wt4_gseReactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

wt7_gseReactome <- gsePathway(WT.geneList[[3]], 
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH", 
                              verbose = FALSE) 
wt7_gseReactome <- setReadable(wt7_gseReactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")


message("+----- Compare clusters gene ontology with BP and Reactome for two approaches------+")
message("+----------wk247vswk0, without split up and dw regulated genes --------------------+")

wtdf1 <- read_excel(outfilesWT[1], sheet = 1); 
wtdf2 <- read_excel(outfilesWT[2], sheet = 1);
wtdf3 <- read_excel(outfilesWT[3], sheet = 1);
wtdf1 <- unique(subset(wtdf1, !is.na(ENTREZID)))
wtdf2 <- unique(subset(wtdf2, !is.na(ENTREZID)))  
wtdf3 <- unique(subset(wtdf3, !is.na(ENTREZID)))  

wtdf1$group <- ifelse(wtdf1$L2FC > 0, "upregulated", "downregulated")
wtdf1$othergroup <- "wk2"
wtdf2$group <- ifelse(wtdf2$L2FC > 0, "upregulated", "downregulated")
wtdf2$othergroup <- "wk4"
wtdf3$group <- ifelse(wtdf3$L2FC > 0, "upregulated", "downregulated")
wtdf3$othergroup <- "wk7"

cwtdf <- rbind(wtdf1, wtdf2, wtdf3)
table(cwtdf$group, cwtdf$othergroup)
#wk2  wk4  wk7
#downregulated  371   95  385
#upregulated   1426  879 1362
## the final analysis may have slightly number difference due to the missingness of geneID,Name or others.

write.csv(cwtdf, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Input_data_compareCluster_wk247vswk0_l2fc06_list_March_2023.csv", row.names=F)

wt_GOAll <- compareCluster(ENTREZID~othergroup, data=cwtdf, fun="enrichGO", OrgDb = org.Mm.eg.db, ont = "BP")
wt_GOAll <- setReadable(wt_GOAll, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
## mouse ID convert to human ID, gProfiler2.
WT.orthosS  <- lapply(WT.orthos, function(x) x[,c(6,7,4)])
wt_orthosDN <- lapply(WT.orthosS, function(x) x[!is.na(x$Human.gene.name),])
wt_orthosGS <- lapply(wt_orthosDN, function(x) bitr(x$Human.gene.name, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db))
## merge the orthos entrezid with L2FC to generate geneList
wt_orthosGL <- list(); wt_genes <- list(); othergroups <- c("wk2", "wk4", "wk7")
for(i in 1:3){
  colnames(wt_orthosDN[[i]]) <- c("ENSEMBL", "SYMBOL", "L2FC")
  wt_orthosGL[[i]] <- merge(wt_orthosDN[[i]], wt_orthosGS[[i]], by = "SYMBOL")
  wt_orthosGL[[i]] <- wt_orthosGL[[i]][order(-wt_orthosGL[[i]]$L2FC),]
  wt_orthosGL[[i]]$group <- ifelse(wt_orthosGL[[i]]$L2FC < 0, "downregulated", "upregulated")
  wt_orthosGL[[i]]$othergroup <- othergroups[i]
}

wt_allGenes <- rbind(wt_orthosGL[[1]], wt_orthosGL[[2]], wt_orthosGL[[3]])

wt_ReactAll <- compareCluster(ENTREZID~othergroup, data=wt_allGenes, fun="enrichPathway", readable=TRUE)

message("+--- Add comparecluster with KEGG pathways, 24/03/2023----------------------+")
wt_KEGGAll <- compareCluster(ENTREZID~othergroup, data=cwtdf, fun="enrichKEGG", organism="mmu")
wt_KEGGAll <- setReadable(wt_KEGGAll, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")



message("+----------wk247vswk0,  split up and dw regulated genes --------------------------+")

wt_GOsplit <- compareCluster(ENTREZID~group+othergroup, data=cwtdf, fun="enrichGO", OrgDb = org.Mm.eg.db, ont = "BP")
wt_GOsplit <- setReadable(wt_GOsplit, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

wt_Reactsplit <- compareCluster(ENTREZID~group+othergroup, data=wt_allGenes, fun="enrichPathway", readable=TRUE)

message("+--- Add comparecluster with KEGG pathways, 24/03/2023----------------------+")
wt_KEGGsplit <- compareCluster(ENTREZID~group+othergroup, data=cwtdf, fun="enrichKEGG", organism="mmu")
wt_KEGGsplit <- setReadable(wt_KEGGsplit, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

save(wt2_reactome, wt4_reactome, wt7_reactome,
     wt2_gseReactome, wt4_gseReactome, wt7_gseReactome,
     wt_GOAll, wt_ReactAll, wt_GOsplit, wt_Reactsplit, 
     wt_KEGGAll, wt_KEGGsplit,
     file="./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_AllReactome_splitBP_Reactome_KEGG_March_2023.RData")

write.xlsx(wt2_reactome, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_AllReactome_March_2023.xlsx",
           row.names=F, append=F, sheetName="wk2vswk0_Reactome")
write.xlsx(wt4_reactome, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_AllReactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk4vswk0_Reactome")
write.xlsx(wt7_reactome, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_AllReactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk7vswk0_Reactome")
write.xlsx(wt2_gseReactome, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_AllReactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk2vswk0_gseReactome")
write.xlsx(wt4_gseReactome, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_AllReactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk4vswk0_gseReactome")
write.xlsx(wt7_gseReactome, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_AllReactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk7vswk0_gseReactome")
####
write.xlsx(wt_GOAll, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_compareCluster_GO_Reactome_March_2023.xlsx",
           row.names=F, append=F, sheetName="wk247vswk0_compareAllGO")
write.xlsx(wt_GOsplit, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_compareCluster_GO_Reactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk247vswk0_compareGOsplit")
write.xlsx(wt_ReactAll, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_compareCluster_GO_Reactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk247vswk0_compareAllReact")
write.xlsx(wt_Reactsplit, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_compareCluster_GO_Reactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk247vswk0_compareReactsplit")
write.xlsx(wt_KEGGAll, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_compareCluster_GO_Reactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk247vswk0_compareAllKEGG")
write.xlsx(wt_KEGGsplit, file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_compareCluster_GO_Reactome_March_2023.xlsx",
           row.names=F, append=T, sheetName="wk247vswk0_compareKEGGsplit")

message("+--------- WT semestic similarity analysis for BP and KEGG--------------------+")
load(gooutfilesWT[1])
allBP.wk2  <- GOSemSim_Data_function(enrichGO_Res,"BP", "wk2vswk0", "allBP", 0.9, Project, 160)
gseBP.wk2  <- GOSemSim_Data_function_GSEBP(enrichGSE_Res, "BP", "wk2vswk0", "gseBP", 0.9, Project, 160)
load(gooutfilesWT[2])
allBP.wk4  <- GOSemSim_Data_function(enrichGO_Res,"BP", "wk4vswk0", "allBP", 0.9, Project, 160)
load(gooutfilesWT[3])
allBP.wk7  <- GOSemSim_Data_function(enrichGO_Res,"BP", "wk7vswk0", "allBP", 0.9, Project, 160)
gseBP.wk7  <- GOSemSim_Data_function_GSEBP(enrichGSE_Res, "BP", "wk7vswk0", "gseBP", 0.9, Project, 160)


message("+--- Check fibrosis relating pathways, including TGF-beta, WNT, hedgehog, Notch, ----+")
message("+--- fibroblast growth factor (FGF) signalling, yes-associated protein1 (YAP)/ ------+")
message("+--- transcriptional coactivator with PDZ binding motif (TAZ) signalling ------------+")
## saved "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx"
## copy from kegg and then save as test.xlsx and then extract gene list.
##FGFR_PI3K <- as.data.frame(read_excel("./Results_Tables_Figures/March_2023/test.xlsx", sheet = 1))
##FGFR_PI3K <- unlist(lapply(FGFR_PI3K[,2], function(x) strsplit(x, split=";")[[1]][1]))

##write.xlsx(as.matrix(TGFbeta_Genes, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
##           append=F, sheetName = "TGFbeta", row.names=F)
##write.xlsx(as.matrix(WNT_Genes, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
##          append=T, sheetName = "WNT", row.names=F)
##write.xlsx(as.matrix(hedgehog_Genes, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
##           append=T, sheetName = "hedgehog", row.names=F)
##write.xlsx(as.matrix(Notch_Genes, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
##          append=T, sheetName = "Notch", row.names=F)
##write.xlsx(as.matrix(FGFR_RAS, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
##           append=T, sheetName = "FGFR_RAS", row.names=F)
##write.xlsx(as.matrix(FGFR_PI3K, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
##          append=T, sheetName = "FGFR_PI3K", row.names=F)
write.xlsx(as.matrix(FGFR_JAKSTAT, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
           append=T, sheetName = "FGFR_JAKSTAT", row.names=F)
write.xlsx(as.matrix(FGFR_MAPK, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
           append=T, sheetName = "FGFR_MAPK1", row.names=F)
write.xlsx(as.matrix(FGFR_PLCGERBB, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
           append=T, sheetName = "FGFR_PLCGERBB", row.names=F)
write.xlsx(as.matrix(FGFR_Hippo, ncol=1), file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Fibrosis_Pathway_GenesList_March_2023.xlsx",
           append=T, sheetName = "FGFR_YAP-TAZ-Hippo", row.names=F)
message("+---Check the overlap for each fibrosis pathways relating with both KEGG and BP for wk4, wk7 kovswt-------+")
## wk4 kovsWT GO and KEGG (enrich)

Fibrosis_GeneList <- list(TGFbeta_Genes, WNT_Genes, hedgehog_Genes, Notch_Genes, FGFR_RAS,
                          FGFR_PI3K, FGFR_JAKSTAT, FGFR_MAPK, FGFR_PLCGERBB, FGFR_Hippo)

## unlist(lapply(Fibrosis_GeneList, length))
## [1]  97 167  57  59 231   359 171 294  84 157
load(gooutfiles1[1])
wk4KO_GOgenes <- lapply(enrichGO_Res$geneID, function(x) strsplit(x, split = "/")[[1]])
overlap_wk4KOGO <- list()
for(i in 1:length(Fibrosis_GeneList)){
  overlap_wk4KOGO[[i]] <- lapply(wk4KO_GOgenes, function(x) intersect(x, Fibrosis_GeneList[[i]]))
}
## wk4/7 formula_KEGGAll
wk47_CompareKEGGAll <- lapply(as.data.frame(formula_KEGGAll)$geneID, function(x) strsplit(x, split = "/")[[1]])
overlap_wk47cpkeggall <- list()
for(i in 1:length(Fibrosis_GeneList)){
  overlap_wk47cpkeggall[[i]] <- lapply(wk47_CompareKEGGAll, function(x) intersect(x, Fibrosis_GeneList[[i]]))
}


message("+-- Only check the KEGG pathways for comparisonclusters all and split up/down-----+")
message("+----                   1) KOvsWT wk4 and wk7                                -----+")
dfKO_All <- formula_KEGGAll %>% 
  mutate(Description = gsub("- Mus musculus \\s*\\([^\\)]+\\)","", Description))
selKEGGAllKO <- c("Viral protein interaction with cytokine and cytokine receptor ", 
                  "Cytokine-cytokine receptor interaction ",
                  "Toll-like receptor signaling pathway ",
                  "Glycine, serine and threonine metabolism ",
                  "Chemokine signaling pathway ",
                  "Biosynthesis of amino acids ",
                  "Cell adhesion molecules ",
                  "IL-17 signaling pathway ",
                  "Diabetic cardiomyopathy ",
                  "Oxidative phosphorylation ",
                  "ECM-receptor interaction ",
                  "Citrate cycle (TCA cycle) ",
                  "Cardiac muscle contraction ",
                  "Hypertrophic cardiomyopathy ",
                  "AGE-RAGE signaling pathway in diabetic complications ",
                  "MAPK signaling pathway ",
                  "NF-kappa B signaling pathway ",
                  "Rap1 signaling pathway ",
                  "B cell receptor signaling pathway ",
                  "Fc gamma R-mediated phagocytosis ",
                  "PI3K-Akt signaling pathway ",
                  "ErbB signaling pathway ",
                  "TNF signaling pathway ",
                  "EGFR tyrosine kinase inhibitor resistance ",
                  "cGMP-PKG signaling pathway ",
                  "Ras signaling pathway ")
pdf("./Results_Tables_Figures/March_2023/Gene_Ontology/KOvsWT_compareKEGG_wk_wk7_all_selPathways_dotplot_28_Mar_2023.pdf",
    width=6, height=12, onefile = T)
dotplot(dfKO_All, showCategory=selKEGGAllKO) +
  theme(axis.text.y = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size = 14, face="bold"),
        axis.title.x = element_text(size = 14, face="bold")) 
dev.off()
##
dfKO_split <- formula_KEGGsplit %>% 
  mutate(Description = gsub("- Mus musculus \\s*\\([^\\)]+\\)","", Description))
pdf("./Results_Tables_Figures/March_2023/Gene_Ontology/KOvsWT_compareKEGG_wk_wk7_split_selPathways_dotplot_28_Mar_2023.pdf",
    width=8, height=12, onefile = T)
dotplot(dfKO_split, showCategory=selKEGGAllKO) +
  theme(axis.text.y = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size = 14, face="bold", angle = 45),
        axis.title.x = element_text(size = 14, face="bold")) 
dev.off()

## Heatmap for selected FGFR pathways for wk7.
Heatmap_FGFR_fn <- function(DEGs_dat, selected_genes, selcolumns, ord.columns, cols.order, pathway, col.labels,
                            name_of_heatmap, condition, gene_number = 50){
  selected_genes <- unique(DEGs_dat[DEGs_dat$external_gene_name %in% selected_genes,])
  selected_genes <- selected_genes[order(selected_genes$log2FoldChange, decreasing = TRUE),] 
  selected_genes <- unique(selected_genes$external_gene_name)
  #if (length(selected_genes) > gene_number )  {
  #  selected_genes <- selected_genes[1:gene_number]
  #}
  rld_mat  <- log2(DEGs_dat[,selcolumns]+1)
  rld_mat1 <- rld_mat - rowMeans(rld_mat)
  rld_mat1$external_gene_name <- DEGs_dat$external_gene_name
  rld_mat1$log2FoldChange <- DEGs_dat$log2FoldChange
  rld_mat2 <- rld_mat1[rld_mat1$external_gene_name%in%selected_genes,]
  rownames(rld_mat2) <- rld_mat2$external_gene_name
  rld_mat2ord <- rld_mat2[match(selected_genes, rld_mat2$external_gene_name), ]   
  rld_mat3 <- subset(rld_mat2ord, select=-c(external_gene_name, log2FoldChange))
  rld_mat3 <- rld_mat3[,ord.columns]
  minval <- round(min(rld_mat3)) - 1
  maxval <- round(max(rld_mat3)) + 1
  f1 = circlize::colorRamp2( c(minval,minval/2,0,maxval/2,maxval), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB") 
  lgd1 = Legend(col_fun = f1, title = "Expression (rld)", at = c(minval,minval/2,0,maxval/2,maxval)  )
  
  L2FC   <- as.matrix(rld_mat2ord$log2FoldChange, ncol = 1)
  pos.num   <- sum(L2FC>0)
  total.len <- dim(L2FC)[1]
  colnames(L2FC) <- "L2FC"
  minfc  <- round(min(L2FC)) - 1
  maxfc  <- round(max(L2FC)) + 1
  fcval  <- max(abs(minfc), abs(maxfc))
  f2 = circlize::colorRamp2( c(minfc,minfc/2,0,maxfc/2,maxfc), c("darkblue", "blue", "grey95", "orange", "darkorange"), space = "RGB") 
  
  ha  = HeatmapAnnotation(Condition = condition,  col = list(Condition = c("WT" = "blue", "KO"= "red")), 
                          annotation_name_side = "left")
  
   
  ht = Heatmap(as.matrix(rld_mat3),  col = f1, name = pathway, show_row_names = T,
               column_title = pathway, cluster_columns = F, row_names_side ="left", cluster_rows = F,
               width = unit(7, "cm"), column_order = cols.order, column_labels = col.labels, 
               show_row_dend = F, row_names_gp = gpar(fontface = "bold.italic"), top_annotation = ha,
               column_names_gp = gpar(fontface = "bold"), column_title_gp = gpar(fontface = "bold"),
               heatmap_legend_param = list(title = "Expression", legend_height = unit(5, "cm"), 
                                           title_position = "leftcenter-rot"),
               row_split = factor(rep(c("A", "B"), c(pos.num, total.len-pos.num))),  row_title = NULL) +
        Heatmap(L2FC,  col = f2, name = "L2FC",  show_row_names = F, width = unit(1, "cm"),
              cluster_columns = F, row_names_side ="right", cluster_rows = F, 
              show_row_dend = F, row_names_gp = gpar(fontface = "bold.italic"), column_names_gp = gpar(fontface = "bold"),
              heatmap_legend_param = list(title = "L2FC", legend_height = unit(5, "cm"), title_position = "leftcenter-rot"),
              row_split = factor(rep(c("A", "B"), c(pos.num, total.len-pos.num))), row_title = NULL)
     
  #pd = packLegend(list = list(lgd1, lgd2, lgd3),direction = "vertical")
  heatmap_height <- nrow(rld_mat3)*0.2+1
  pdf(paste("./Results_Tables_Figures/March_2023/Gene_Ontology/Fig_xx_", "ComplexHeatmap_",  name_of_heatmap, ".pdf", sep=""), onefile=TRUE, width=9, height=heatmap_height) 
  par(bg=NA)
  draw(ht, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", 
       gap = unit(1, "cm"), legend_grouping = "original")
  dev.off()
  
}

datwk7_counts <- read.csv("./Results_Tables_Figures/March_2023/CAD_jcg23_0001-Resall_lfcShrink_7w_summaryTable_March_2023.csv", header = T)
c <- datwk7_counts
condition <- rep(c("WT", "KO"), each=3)
formula_KEGGAll <- dfKO_All
selected_genes <- unlist(strsplit(formula_KEGGAll[formula_KEGGAll@compareClusterResult$Description=="PI3K-Akt signaling pathway ", "geneID"], split = "/"))
selcolumns <- 14:19
ord.columns <- c(4:6, 1:3)
cols.order <- c(1:6)
pathway <- "PI3K-Akt signaling pathway"
col.labels <- c("WT_R1", "WT_R2", "WT_R3", "KO_R1", "KO_R2", "KO_R3")
name_of_heatmap <- "PI3K_wk7_KOvsWT"
PI3K <- Heatmap_FGFR_fn(DEGs_dat, selected_genes, selcolumns, ord.columns, cols.order, pathway, col.labels,
                                    name_of_heatmap, condition, gene_number = length(selected_genes))
selected_genes <- unlist(strsplit(formula_KEGGAll[formula_KEGGAll@compareClusterResult$Description=="Ras signaling pathway ", "geneID"], split = "/"))
pathway <- "Ras signaling pathway "
name_of_heatmap <- "RAS_wk7_KOvsWT"
RAS <- Heatmap_FGFR_fn(DEGs_dat, selected_genes, selcolumns, ord.columns, cols.order, pathway, col.labels,
                        name_of_heatmap, condition, gene_number = length(selected_genes))
selected_genes <- unlist(strsplit(formula_KEGGAll[formula_KEGGAll@compareClusterResult$Description=="MAPK signaling pathway ", "geneID"], split = "/"))
pathway <- "MAPK signaling pathway "
name_of_heatmap <- "MAPK_wk7_KOvsWT"
MAPK <- Heatmap_FGFR_fn(DEGs_dat, selected_genes, selcolumns, ord.columns, cols.order, pathway, col.labels,
                       name_of_heatmap,condition, gene_number = length(selected_genes))
selected_genes <- unlist(strsplit(formula_KEGGAll[formula_KEGGAll@compareClusterResult$Description=="ErbB signaling pathway ", "geneID"], split = "/"))
pathway <- "ErbB signaling pathway "
name_of_heatmap <- "ErbB_wk7_KOvsWT"
ERBB <- Heatmap_FGFR_fn(DEGs_dat, selected_genes, selcolumns, ord.columns,cols.order, pathway, col.labels,
                        name_of_heatmap, condition,gene_number = length(selected_genes))

message("+---- up/dw split KEGG pathways-------------+")
dfKO_split <- formula_KEGGsplit %>% 
  mutate(Description = gsub("- Mus musculus \\s*\\([^\\)]+\\)","", Description))
dotplot(dfKO_split, showCategory =70)
selKEGGsplitKO <- c(selKEGGAllKO, "TGF-beta signaling pathway ", "p53 signaling pathway " ,
                    "NOD-like receptor signaling pathway ","C-type lectin receptor signaling pathway ")
dfKO_splitdat <- as.data.frame(dfKO_split)
selHeatKEGGsplit <- dfKO_splitdat[dfKO_splitdat$Description%in%selKEGGsplitKO,]                  
selHeatPltKEGGsplit <- selHeatKEGGsplit[,c("Cluster", "Description", "GeneRatio", "p.adjust")]


selHeatMatplt <- matrix(0, ncol = 4, nrow=27)
for(i in 1:length(unique(selHeatPltKEGGsplit$Description))){
  selHeatMatplt[i,1] <- ifelse(selHeatPltKEGGsplit$Cluster=="downregulated.wk4" & selHeatPltKEGGsplit$Description == unique(selHeatPltKEGGsplit$Description)[i], selHeatPltKEGGsplit$p.adjust, "NA")
  selHeatMatplt[i,2] <- ifelse(selHeatPltKEGGsplit$Cluster=="upregulated.wk4" & selHeatPltKEGGsplit$Description == unique(selHeatPltKEGGsplit$Description)[i], selHeatPltKEGGsplit$p.adjust, "NA")
  selHeatMatplt[i,3] <- ifelse(selHeatPltKEGGsplit$Cluster=="downregulated.wk7" & selHeatPltKEGGsplit$Description == unique(selHeatPltKEGGsplit$Description)[i], selHeatPltKEGGsplit$p.adjust, "NA")
  selHeatMatplt[i,4] <- ifelse(selHeatPltKEGGsplit$Cluster=="upregulated.wk7" & selHeatPltKEGGsplit$Description == unique(selHeatPltKEGGsplit$Description)[i], selHeatPltKEGGsplit$p.adjust, "NA")
  
  }

message("+--- FINISH 15/03/2023---------------------+")

message("+---- DESeq2 analysis based on KO type across the time-------------------------+")

resultsNames(dds_KO)
resKO.2wvs0w      <- lfcShrink(dds_KO, coef="Time_2w_vs_0w", type="apeglm")
resKO.2wvs0w.sig  <- as.data.frame(subset(resKO.2wvs0w, padj <= significance & abs(log2FoldChange) >= l2fc))
print(nrow(resKO.2wvs0w.sig)); print(table(sign(resKO.2wvs0w.sig$log2FoldChange)))
resKO.2wvs0w.sig$ensembl_gene_id <- rownames(resKO.2wvs0w.sig)
resKO.2wvs0w.sigmer <- merge(resKO.2wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
print(nrow(resKO.2wvs0w.sig)); print(table(sign(resKO.2wvs0w.sig$log2FoldChange)))
resKO.2wvs0w.sig$ensembl_gene_id <- rownames(resKO.2wvs0w.sig)
resKO.2wvs0w.sigmer <- merge(resKO.2wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
write.xlsx(resKO.2wvs0w.sigmer, file = "./Results_Tables_Figures/March_2023/KO_timemodel_wk247vs0_sigDEGs_l2fc1_Mar_2023.xlsx",
           row.names=F, append = F, sheetName= "KO_2w_vs_0w_sig")
RESKO.files <- paste0("./Results_Tables_Figures/March_2023/KO_timemodel_wk", c(2,4, 7), "vs0_allDEGs_Mar_2023.csv")
resKO.2wvs0w.merEC <- merECRes_function_pair(resKO.2wvs0w, ensEMBL2id, dds_KO, RESKO.files[1])

##
resKO.4wvs0w      <- lfcShrink(dds_KO, coef="Time_4w_vs_0w", type="apeglm")
resKO.4wvs0w.sig  <- as.data.frame(subset(resKO.4wvs0w, padj <= significance & abs(log2FoldChange) >= l2fc))
print(nrow(resKO.4wvs0w.sig)); print(table(sign(resKO.4wvs0w.sig$log2FoldChange)))
resKO.4wvs0w.sig$ensembl_gene_id <- rownames(resKO.4wvs0w.sig)
resKO.4wvs0w.sigmer <- merge(resKO.4wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
print(nrow(resKO.4wvs0w.sig)); print(table(sign(resKO.4wvs0w.sig$log2FoldChange)))
resKO.4wvs0w.sig$ensembl_gene_id <- rownames(resKO.4wvs0w.sig)
resKO.4wvs0w.sigmer <- merge(resKO.4wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
write.xlsx(resKO.4wvs0w.sigmer, file = "./Results_Tables_Figures/March_2023/KO_timemodel_wk247vs0_sigDEGs_l2fc1_Mar_2023.xlsx",
           row.names=F, append = T, sheetName= "KO_4w_vs_0w_sig")
resKO.4wvs0w.merEC <- merECRes_function_pair(resKO.4wvs0w, ensEMBL2id, dds_KO, RESKO.files[2])

##
resKO.7wvs0w      <- lfcShrink(dds_KO, coef="Time_7w_vs_0w", type="apeglm")
resKO.7wvs0w.sig  <- as.data.frame(subset(resKO.7wvs0w, padj <= significance & abs(log2FoldChange) >= l2fc))
print(nrow(resKO.7wvs0w.sig)); print(table(sign(resKO.7wvs0w.sig$log2FoldChange)))
resKO.7wvs0w.sig$ensembl_gene_id <- rownames(resKO.7wvs0w.sig)
resKO.7wvs0w.sigmer <- merge(resKO.7wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
print(nrow(resKO.7wvs0w.sig)); print(table(sign(resKO.7wvs0w.sig$log2FoldChange)))
resKO.7wvs0w.sig$ensembl_gene_id <- rownames(resKO.7wvs0w.sig)
resKO.7wvs0w.sigmer <- merge(resKO.7wvs0w.sig, ensEMBL2id, by ="ensembl_gene_id")
write.xlsx(resKO.7wvs0w.sigmer, file = "./Results_Tables_Figures/March_2023/KO_timemodel_wk247vs0_sigDEGs_l2fc1_Mar_2023.xlsx",
           row.names=F, append = T, sheetName= "KO_7w_vs_0w_sig")
resKO.7wvs0w.merEC <- merECRes_function_pair(resKO.7wvs0w, ensEMBL2id, dds_KO, RESKO.files[3])


message("+----------------------------------------------------------------------------------------------------------+")
message("+--- Generate a plot with the number of sigDEGs across time for WT, KO plus KOvsWT            -------------+")
message("+----------------------------------------------------------------------------------------------------------+")

RES.files1  <- paste0("./Results_Tables_Figures/March_2023/",Project, "-Resall_lfcShrink_",Time, "_summaryTable_March_2023.csv")
RESWT.files <- paste0("./Results_Tables_Figures/March_2023/WT_timemodel_wk", c(2,4, 7), "vs0_allDEGs_Mar_2023.csv")
RESKO.files <- paste0("./Results_Tables_Figures/March_2023/KO_timemodel_wk", c(2,4, 7), "vs0_allDEGs_Mar_2023.csv")

DEGs_Table_fn <- function(Resfile, significance, l2fc){
  resdat <- read.csv(Resfile, header = T)
  subdim <- nrow(subset(resdat, padj <= significance & abs(log2FoldChange) >= l2fc))
  subpos <- nrow(subset(resdat, padj <= significance & log2FoldChange >= l2fc))
  subneg <- nrow(subset(resdat, padj <= significance & log2FoldChange <= -l2fc))
  subs <- c(subdim, subpos, subneg)
  subs
}
DEcount_Table_KOvsWT_1 <- lapply(RES.files1, function(x) DEGs_Table_fn(x, 0.05, 1))
DEcount_Table_KOvsWT_06 <- lapply(RES.files1, function(x) DEGs_Table_fn(x, 0.05, 0.6))
DEcount_Table_WT_1 <- lapply(RESWT.files, function(x) DEGs_Table_fn(x, 0.05, 1))
DEcount_Table_WT_06 <- lapply(RESWT.files, function(x) DEGs_Table_fn(x, 0.05, 0.6))
DEcount_Table_KO_1 <- lapply(RESKO.files, function(x) DEGs_Table_fn(x, 0.05, 1))
DEcount_Table_KO_06 <- lapply(RESKO.files, function(x) DEGs_Table_fn(x, 0.05, 0.6))

DE_Count_Table <- melt( t( data.frame("wk0_KOvsWT"=DEcount_Table_KOvsWT_1[[1]][1], 
                                      "wk2_KOvsWT"=DEcount_Table_KOvsWT_1[[2]][1], 
                                      "wk4_KOvsWT"=DEcount_Table_KOvsWT_1[[3]][1],
                                      "wk7_KOvsWT"=DEcount_Table_KOvsWT_1[[4]][1],
                                      "WT_wk2vswk0"=DEcount_Table_WT_1[[1]][1],
                                      "WT_wk4vswk0"=DEcount_Table_WT_1[[2]][1],
                                      "WT_wk7vswk0"=DEcount_Table_WT_1[[3]][1],
                                      "KO_wk2vswk0"=DEcount_Table_KO_1[[1]][1],
                                      "KO_wk4vswk0"=DEcount_Table_KO_1[[2]][1],
                                      "KO_wk7vswk0"=DEcount_Table_KO_1[[3]][1]) ) )

DE_Count_Table$Var2 <- rep(c("KOvsWT", "WT", "KO"), c(4,3,3))
DE_Count_Table$Time <- c("wk0", "wk2", "wk4", "wk7", "wk2", "wk4", "wk7","wk2", "wk4", "wk7")
DE_Count_Table$Week <- factor(gsub("wk", "", DE_Count_Table$Time),levels=c("0","2", "4","7"))
colnames(DE_Count_Table) <- c("Sample", "Condition", "Counts",  "Time", "Week")

## Refined the color scheme consistent with PCA plot.
#pdf(paste0("./Results_Tables_Figures/March_2023/", Project, "_Fig.Comparing_Numbers_DEGs_with_time_April_2023.pdf"),width=5,height=7)
pdf(paste0("./Results_Tables_Figures/Paper_v1_April_2023/SuppFigxx-Comparing_Numbers_DEGs_with_time_April_2023.pdf"),width=5,height=4)
par(bg=NA)
ggplot(data=DE_Count_Table, aes(x=Week, y=Counts, group=Condition, colour=Condition)) +
  geom_point(size= 3, alpha=0.75) +
  geom_line(linewidth=1.2, alpha=0.75) +
  scale_colour_manual(name="Treatment", values=c("KOvsWT"="purple", "WT"="royalblue",  "KO"="darkorange")) +
  theme(axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=14, face="bold"),
        axis.title.y=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=10, face= "bold")) + 
  ylab("Number of DEGs") +
  scale_x_discrete(name="Week", limits=c("0", "2","4","7"))
dev.off()



message("+----------------------------------------------------------------------------------------------------------+")
message("+--- Select the key markers relating to the inflamsome,fibrosis, metabolic, serine-1 pathways--------------+")
message("+----------------------------------------------------------------------------------------------------------+")

message("+-----------------Retrieve ensEMBL annotations-----------------+")

load("/Users/xz289/Documents/Minoru/Project_mt709_0001/Data/References_Data/GRCm39/Ensembl_GRCm39_entrezID_ensembID_exterName.RData")
nrow(ensEMBL2id)
load("/Users/xz289/Documents/Minoru/Project_mt709_0001/Data/References_Data/GRCm39/Ensembl_mmusculus_ID_Chr_GOID_GOName_GRCm39.RData")
nrow(ensEMBL2id.GOInfo)
ensEMBL2id_gomer <- merge(ensEMBL2id, ensEMBL2id.GOInfo, by = "ensembl_gene_id")[,-7]
colnames(ensEMBL2id_gomer)[5] <-"chromosome"
ensEMBL2id_gomer$description  <- gsub("\\[.*?\\]","",as.character(ensEMBL2id_gomer$description))
save(ensEMBL2id_gomer, file="/Users/xz289/Documents/Minoru/Project_mt709_0001/Data/References_Data/GRCm39/Ensembl_GRCm39_entrezID_ensembID_exterName_GOID_GOName.RData")


message("+-----------------Select the corresponding markers with our interested pathways-----------------+")

# inflammasome genes:
ensEMBL2id_go_Inflammasome     <- subset(ensEMBL2id_gomer, grepl("inflammasome", ensEMBL2id_gomer$name_1006))
Inflammasome_GO_terms          <- unique(ensEMBL2id_go_Inflammasome$name_1006)
Inflammasome_GO_genes_ensembl  <- unique(ensEMBL2id_go_Inflammasome$ensembl_gene_id)
Inflammasome_GO_genes_external <- unique(ensEMBL2id_go_Inflammasome$external_gene_name)
length(Inflammasome_GO_genes_ensembl)

# inflammatory genes:
ensEMBL2id_go_Inflammatory     <- subset(ensEMBL2id_gomer, grepl("inflammatory", ensEMBL2id_gomer$name_1006))
Inflammatory_GO_terms          <- unique(ensEMBL2id_go_Inflammatory$name_1006)
Inflammatory_GO_genes_external <- unique(ensEMBL2id_go_Inflammatory$external_gene_name)
Inflammatory_GO_genes_ensembl  <- unique(ensEMBL2id_go_Inflammatory$ensembl_gene_id)
length(Inflammatory_GO_genes_ensembl)

Inflammatory_GO_genes_external[Inflammatory_GO_genes_external%in%Inflammasome_GO_genes_external] ## 22 common genes
## "Ddx3x"  "Naip1"  "Mefv"   "Gsdmd"  "Casp12" "Casp1"  "Pycard" "Nlrp3"  "Casp4"  "Aim2"   "Nlrp6"  "Tlr4"   "Nlrc4" 
## "Nlrc3"  "Tlr6"   "Nlrp9b" "Nlrp1a" "Nlrp1b" "Naip5"  "Naip6"  "Naip2"  "Gbp5"  

# Fibrosis genes:https://www.informatics.jax.org/vocab/gene_ontology/ fibrosis (wound healing)
ensEMBL2id_go_Fibrosis         <- subset(ensEMBL2id_gomer, (ensEMBL2id_gomer$go_id == "GO:0002248" | 
                                                            ensEMBL2id_gomer$go_id == "GO:1904596" | 
                                                            ensEMBL2id_gomer$go_id == "GO:0097709" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904597" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904598" ) )
fibrosis_GO_terms              <- unique(ensEMBL2id_go_Fibrosis$name_1006)
fibrosis_GO_genes_external     <- unique(ensEMBL2id_go_Fibrosis$external_gene_name)
fibrosis_GO_genes_ensembl      <- unique(ensEMBL2id_go_Fibrosis$ensembl_gene_id)
length(fibrosis_GO_genes_ensembl)
# myofibroblast
ensEMBL2id_go_myofibroblast    <- subset(ensEMBL2id_gomer, (ensEMBL2id_gomer$go_id == "GO:0036446" | 
                                                            ensEMBL2id_gomer$go_id == "GO:1904516" | 
                                                            ensEMBL2id_gomer$go_id == "GO:1990764" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904329" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904330" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904521" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904522" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904760" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904761" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904762" |
                                                            ensEMBL2id_gomer$go_id == "GO:1904328"  ) )
myofibroblast_GO_terms              <- unique(ensEMBL2id_go_myofibroblast$name_1006)
myofibroblast_GO_genes_external     <- unique(ensEMBL2id_go_myofibroblast$external_gene_name)
myofibroblast_GO_genes_ensembl      <- unique(ensEMBL2id_go_myofibroblast$ensembl_gene_id)
length(myofibroblast_GO_genes_ensembl)


Fibrosis_myofibroblast <- unique(c("Tgfb1","Acta2","Smad4","Smad2","Mir669i","Tmem9b","Cask","Smad3","Gm23336",
"Tgfbr2","Pdilt","Ctbp2","Tgfbr1","Stim1","Bccip","Dhx32", myofibroblast_GO_genes_external))

# Serine-1 (metabolic dehydrogenase)
ensEMBL2id_go_serine1     <- subset(ensEMBL2id_gomer, grepl("L-serine", ensEMBL2id_gomer$name_1006))
serine1_GO_terms          <- unique(ensEMBL2id_go_serine1$name_1006)
serine1_GO_genes_external <- unique(ensEMBL2id_go_serine1$external_gene_name)
serine1_GO_genes_ensembl  <- unique(ensEMBL2id_go_serine1$ensembl_gene_id)
length(serine1_GO_genes_ensembl)

# ATF4
ensEMBL2id_go_atf4     <- subset(ensEMBL2id_gomer, grepl("ATF4", ensEMBL2id_gomer$name_1006))
atf4_GO_terms          <- unique(ensEMBL2id_go_atf4$name_1006)
atf4_GO_genes_external <- unique(ensEMBL2id_go_atf4$external_gene_name)
atf4_GO_genes_ensembl  <- unique(ensEMBL2id_go_atf4$ensembl_gene_id)
length(atf4_GO_genes_ensembl)

# metabolic
ensEMBL2id_go_metabolic <- subset(ensEMBL2id_gomer, grepl("metabolic", ensEMBL2id_gomer$name_1006))
metabolic_GO_terms          <- unique(ensEMBL2id_go_metabolic$name_1006)
metabolic_GO_genes_external <- unique(ensEMBL2id_go_metabolic$external_gene_name)
metabolic_GO_genes_ensembl  <- unique(ensEMBL2id_go_metabolic$ensembl_gene_id)
length(metabolic_GO_genes_ensembl)

# p53 
ensEMBL2id_go_p53 <- subset(ensEMBL2id_gomer, grepl("p53", ensEMBL2id_gomer$name_1006))
p53_GO_terms          <- unique(ensEMBL2id_go_p53$name_1006)
p53_GO_genes_external <- unique(ensEMBL2id_go_p53$external_gene_name)
p53_GO_genes_ensembl  <- unique(ensEMBL2id_go_p53$ensembl_gene_id)
length(p53_GO_genes_ensembl)


message("+--- The above category selected markers relating to DEGs list and pathways. mainly for wk4 and wk7 KOvsWT.---+")

Fibrosis_Keywords  <- c("wound healing", "PI3K-Akt", "Ras", "MAPK", "ErbB")
Metabolic_Keywords <- c("amino acid metabolic", "fatty acid", "collagen", "Glycine, serine", "Arginine", "Cytokine-cytokine")
Inflam_Keywords    <- c("inflammatory", "inflammasome")
SerineOne_Keywords <- c("L-serine", "Biosynthesis of amino acids")
ATF4_Keywords      <- c("PERK", "ATF4", "EIF2AK1")
stress_Keywords    <- c("stress") # ISR includes ATF4, PERK, CHOP...
CHOP_Keywords      <- c("apoptosis", "cell death", "granzyme A", "DNA damage")
EM_Keywords        <- "extracellular"
Cellcycle_Keywords <- c("Cell cycle", "TCA")
P53_Keywords       <- c("p53", "TP53")
# just loaded the comparison GO/KEGG/Reactome pathways files.
load("./Results_Tables_Figures/March_2023/Gene_Ontology/CompareCluster_wk47_KOvsWT_all_splitupdw_April_2023.RData")

pathway_KeywordsList <- list(Fibrosis_Keywords, Metabolic_Keywords, Inflam_Keywords, SerineOne_Keywords,
                             ATF4_Keywords, stress_Keywords, CHOP_Keywords, EM_Keywords, Cellcycle_Keywords,
                             P53_Keywords)
GOsearch_function <- function(GOinput, keywordsList){
  GOdat  <- as.data.frame(GOinput)
  Outdat <- lapply(keywordsList, function(x) GOdat[grep(x, GOdat$Description), ])
  Outdat <- Outdat[sapply(Outdat, nrow) != 0]
  Outdat <- ldply(Outdat, data.frame) 
  Outdat
}
SelPathNames <- c("Fibrosis", "Metabolic", "Inflam", "SerineOne", "ATF4", "stress","CHOP", "EM", "Cellcycle", "P53")
GOKO_outputNames   <- list(paste0(SelPathNames, "_GOKO"))
KEGGKO_outputNames <- list(paste0(SelPathNames, "_KEGGKO"))
ReactKO_outputNames <- list(paste0(SelPathNames, "_ReactKO"))
for(i in 1:length(pathway_KeywordsList)){
  GOKO_outputNames[[i]]    <- GOsearch_function(formula_GOAll, pathway_KeywordsList[[i]])
  KEGGKO_outputNames[[i]]  <- GOsearch_function(formula_KEGGAll, pathway_KeywordsList[[i]])
  ReactKO_outputNames[[i]] <- GOsearch_function(formula_ReactAll, pathway_KeywordsList[[i]])
}

write.xlsx(GOKO_outputNames[[1]], file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Selected_wk4_wk7_KOvsWT_Pathways_summary_18_April_2023.xlsx", 
           append=F, sheetName = paste0(SelPathNames[1], "_GOKO"))
for(i in 1:10){
  print(i)
  if(dim(KEGGKO_outputNames[[i]])[1]!=0){
    write.xlsx(KEGGKO_outputNames[[i]], file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Selected_wk4_wk7_KOvsWT_Pathways_summary_18_April_2023.xlsx", 
             append=T, sheetName = paste0(SelPathNames[i], "_KEGGKO"))}
  if(dim(ReactKO_outputNames[[i]])[1]!=0){
    write.xlsx(ReactKO_outputNames[[i]], file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Selected_wk4_wk7_KOvsWT_Pathways_summary_18_April_2023.xlsx", 
             append=T, sheetName = paste0(SelPathNames[i], "_ReactKO"))}
  if(dim(GOKO_outputNames[[i+1]])[1]!=0){
    write.xlsx(GOKO_outputNames[[i+1]], file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Selected_wk4_wk7_KOvsWT_Pathways_summary_18_April_2023.xlsx", 
            append=T, sheetName = paste0(SelPathNames[i+1], "_GOKO"))}
}
message("+--- after selection and to generate heatmap for wk4/7 KOvsWT pathways-------------------+")

selKOpathways <- readxl::read_excel("./Results_Tables_Figures/March_2023/Gene_Ontology/Selpathway_Heatmap_list_April_2023.xlsx")
selKOpathways$adjust.pval <- -log10(selKOpathways$p.adjust)
library(tidyverse)

#Pipe and merge
dataset2 <- selKOpathways %>%
  group_by(Description) %>%
  mutate(id_rows = othergroup) %>%
  pivot_wider( 
    id_cols = c(Description, 'GO Category', ParentTerm),
    names_from = id_rows,
    values_from = c(adjust.pval)
    
  ) %>% 
  ungroup()
dataset2 <- as.data.frame(dataset2)
dataset2$ParentTerm <- gsub("CellCycle", "TCA cycle", dataset2$ParentTerm)
dataset2$ParentTerm <- gsub("Inflammatory", "Inflam", dataset2$ParentTerm)
rownames(dataset2) <- dataset2$Description
library(RColorBrewer)
cols <- brewer.pal(6, "Greens")
pal  <- colorRampPalette(cols)
ht_cols = pal(20)
rownames(dataset2) <- c("wound healing", "wound healing involved in \ninflammatory response",
                        "vascular wound healing", "PI3K-Akt signaling pathway",
                        "Ras signaling pathway", "MAPK signaling pathway",
                        "ErbB signaling pathway", "glutamine family amino \nacid metabolic process",
                        "alpha-amino acid metabolic process", "cellular amino acid metabolic process",
                        "serine family amino acid \nmetabolic process","aspartate family amino \nacid metabolic process",
                        "collagen metabolic process", "collagen fibril organization", "collagen biosynthetic process",
                        "Glycine, serine and threonine metabolism","Arginine and proline metabolism", 
                        "Cytokine-cytokine receptor interaction", "regulation of inflammatory response",
                        "acute inflammatory response", "Cell recruitment \n(pro-inflammatory response)",
                        "Inflammasomes", "L-serine metabolic process", "Biosynthesis of amino acids",
                        "PERK-mediated unfolded protein response", "PERK regulates gene expression",
                        "ATF4 activates genes in response \nto endoplasmic reticulum  stress",
                        "Response of EIF2AK1 (HRI) to heme deficiency","regulation of translation \nin response to stress", 
                        "cellular response to chemical stress", "response to oxidative stress", "response to endoplasmic reticulum stress",
                        "intrinsic apoptotic signaling pathway in \nresponse to endoplasmic reticulum stress",
                        "positive regulation of stress-activated \nMAPK cascade",
                        "intrinsic apoptotic signaling pathway \nin response to DNA damage",
                        "response to extracellular stimulus", "extracellular matrix organization",
                        "cellular response to extracellular stimulus", "The citric acid (TCA) cycle and \nrespiratory electron transport",
                        "Pyruvate metabolism and \nCitric Acid (TCA) cycle", "Citric acid cycle (TCA cycle)",
                        "positive regulation of signal \ntransduction by p53 class mediator", "signal transduction by p53 class mediator",
                        "intrinsic apoptotic signaling pathway in response \nto DNA damage by p53 class mediator",
                        "intrinsic apoptotic signaling pathway by p53 class mediator","regulation of signal transduction by p53 class mediator",
                        "DNA damage response, signal transduction by p53 class mediator")


pdf("./Results_Tables_Figures/March_2023/Gene_Ontology/test_Heatmap_selKOpathways_April_2023.pdf", height= 12, width=12)
KOheat1 <- ComplexHeatmap:: Heatmap(as.matrix(dataset2[,c("wk4","wk7")]),  
                         col = ht_cols, name = "Pathways",  
                         column_title = "", 
                         show_row_names = TRUE, show_column_names = TRUE, 
                         heatmap_legend_param = list(title = "-log10(padj)", 
                                                     legend_height = unit(3, "cm"), 
                                                     title_position = "topleft"), 
                         cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left", 
                         column_names_rot = 0, column_names_side = "bottom", 
                         column_title_rot = 0, column_names_centered = TRUE,
                         column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                         row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                         column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                         row_split = dataset2$ParentTerm, row_gap = unit(5, "mm"), 
                         na_col = "darkgray", row_title = NULL,
                         width = ncol(dataset2)*unit(4, "mm"), 
                         height = nrow(dataset2)*unit(7, "mm"))
dev.off()

message("+--- selected the same pathways from wk2/4,7 vs wk0 for WT geneontology analysis ---+")

load("./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_AllReactome_splitBP_Reactome_KEGG_March_2023.RData")
selGOID   <- unique(selKOpathways$ID)
selDat1   <- unique(selKOpathways[c("Description", "GO Category", "ParentTerm")]) ## 41
WT_GOsels <- as.data.frame(wt_GOAll)[as.data.frame(wt_GOAll)$ID%in%selGOID, ] # 41
WT_KEGGsels  <- as.data.frame(wt_KEGGAll)[as.data.frame(wt_KEGGAll)$ID%in%selGOID, ] # 7
WT_Reactsels <- as.data.frame(wt_ReactAll)[as.data.frame(wt_ReactAll)$ID%in%selGOID, ] # 2
WT_GOmer1 <- rbind(WT_GOsels, WT_KEGGsels,WT_Reactsels) ## unique 23
WT_GOmer1$adjust.pval <- -log10(WT_GOmer1$p.adjust)
WT_GOmer2 <- merge(selDat1, WT_GOmer1, by = "Description", all.x=T)
dataset_wt <-WT_GOmer2 %>%
  group_by(Description) %>%
  mutate(id_rows = othergroup) %>%
  pivot_wider( 
    id_cols = c(Description, 'GO Category', ParentTerm),
    names_from = id_rows,
    values_from = c(adjust.pval)
 ) %>% 
  ungroup()
dataset_wt <- dataset_wt[, -5]
dataset_wt <- dataset_wt[,c(1,2,3,5,6,4)]

order.levels <- factor(dataset_wt$Description, levels = as.factor(selDat1$Description))
dataset_wt <- as.data.frame(dataset_wt)[order(order.levels),]
rownames(dataset_wt) <- dataset_wt$Description
## 
cols1 <- brewer.pal(6, "Purples")
pal1  <- colorRampPalette(cols1)
ht_cols1 = pal1(20)

row_ha = rowAnnotation(Class= dataset_wt$ParentTerm)
WTheat1 <- ComplexHeatmap:: Heatmap(as.matrix(dataset_wt[,c("wk2","wk4","wk7")]),  
                                    col = ht_cols1, name = "WTPathways",  
                                    row_title = "", column_title = "", 
                                    show_row_names = TRUE, show_column_names = TRUE, 
                                    heatmap_legend_param = list(title = "-log10(padj)", 
                                                                legend_height = unit(3, "cm"), 
                                                               title_position = "topleft"), 
                                    #show_heatmap_legend = FALSE,
                                    cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left", 
                                    column_names_rot = 0, column_names_side = "bottom", 
                                    column_title_rot = 0, column_names_centered = TRUE,
                                    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                    row_split = dataset2$ParentTerm, row_gap = unit(5, "mm"), row_title_side ="right",
                                    row_title_gp = gpar(fontsize = 12, fontface = "bold"), na_col = "darkgray",
                                    width = ncol(dataset2)*unit(4, "mm"), 
                                    height = nrow(dataset2)*unit(7, "mm"))

pdf("./Results_Tables_Figures/March_2023/Gene_Ontology/Figxx-Pathways_Heatmap_KOvsWT_WT_19_April_2023_addp53.pdf",
    height= 14, width=12, onefile=T)
KOheat1+WTheat1+row_ha
dev.off()
# need to make the legend with one not two. 

message("+---- Generate heatmap with selfdefined groups pathways using the padj mean across above-----+")
gowt_dat <- dataset_wt; gowt_dat[is.na(gowt_dat)] <- 0
goko_dat <- dataset2; goko_dat[is.na(goko_dat)] <- 0

gowt_means <- aggregate(. ~ ParentTerm, data = gowt_dat[,c(3:6)], FUN = mean, na.rm=T)
goko_means <- aggregate(. ~ ParentTerm, data = goko_dat[,c(3:5)], FUN = mean, na.rm=T)
gowt_means[2,1] <- "TCA cycle"; gowt_means[5,1] <- "Inflam"
gowt_means <- gowt_means[order(gowt_means[,1]), ];
goko_means <- goko_means[order(goko_means[,1]), ];
groupNames <- c("ATF4 activates in response to ER stress and ISR",
                "Extracellualr matrix organisation",
                "Fibrosis",
                "Inflamsome/inflammation relating", 
                "Metabolic relating",
                "P53 signalling pathway",
                "Serine biosynthesis and one carbon cycle",
                "Response to stress",
                "Pyruvate metabolism and TCA cycle")
rownames(gowt_means) <- groupNames;
rownames(goko_means) <- groupNames
gowt_means <- gowt_means[,-1]; gowt_means[gowt_means==0] <- NA
goko_means <- goko_means[,-1]; goko_means[goko_means==0] <- NA

cols2 <- brewer.pal(4, "Greens")
pal2  <- colorRampPalette(cols2)
ht_cols2 = pal2(20)

gowt_hplt <- ComplexHeatmap:: Heatmap(as.matrix(gowt_means),  
                         col = ht_cols2,
                         name = "WTPathways",  
                         row_title = "", column_title = "WT2/4/7 \nvs \nWT0", 
                         show_row_names = TRUE, show_column_names = TRUE,
                         column_labels = c(2,4,7),
                         heatmap_legend_param = list(title = "-log10(padj)", 
                                                     legend_height = unit(3, "cm"), 
                                                     title_position = "leftcenter-rot"), 
                         #show_heatmap_legend = FALSE,
                         cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left", 
                         column_names_rot = 0, column_names_side = "bottom", na_col = "darkgray",
                         column_title_rot = 0, column_names_centered = TRUE,
                         column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                         row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                         column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                         width = ncol(gowt_means)*unit(7, "mm"), 
                         height = nrow(gowt_means)*unit(7, "mm"))

goko_hplt <- ComplexHeatmap:: Heatmap(as.matrix(goko_means[,c("wk4","wk7")]),  
                                      col = ht_cols2, name = "KOvsWTPathways",  
                                      row_title = "", column_title = "KO vs WT", 
                                      show_row_names = TRUE, show_column_names = TRUE,
                                      column_labels = c(4,7), 
                                      heatmap_legend_param = list(title = "-log10(padj)", 
                                                                  legend_height = unit(3, "cm"), 
                                                                  title_position = "leftcenter-rot"), 
                                      show_heatmap_legend = FALSE,
                                      cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left", 
                                      column_names_rot = 0, column_names_side = "bottom", 
                                      column_title_rot = 0, column_names_centered = TRUE,
                                      column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                      row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                      column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                      width = ncol(gowt_means)*unit(5, "mm"), 
                                      height = nrow(gowt_means)*unit(7, "mm"))
pdf("./Results_Tables_Figures/Paper_v1_April_2023/Figxx-Pathways_Heatmap_KOvsWT_WT_avgP_25_April_2023.pdf",
    height= 6, width=10, onefile=T)
goko_hplt + gowt_hplt
dev.off()

message("+-----------------Bubble plot for the above markers   -----------------+")

makeCustomGeneSetBubbleplots <- function(Base.dir, GENELIST, GENENAME, DDS, ensEMBL2id, Fig.Height, Fig.Width, 
                                         Fig.ncol, TITLE, textsize, axissize, titlesize) {
  
  count            <- 0
  BigTable         <- data.frame(count=double(), condition=character(), gene=character(), genename=character())
  BigTable.summary <- data.frame(condition=character(), mean=double(), gene=character(), genename=character())
  
  for(gene in GENELIST)  
  {
    message(gene)
   # colData(DDS)$condition <- paste0(colData(DDS)$Condition, "_", colData(DDS)$Time)
    tmp.TC          <- plotCounts(DDS, gene, normalised=TRUE, intgroup = c("Condition", "Time"), returnData = TRUE)
    tmp.TC$gene     <- gene
    tmp.TC$genename <- subset(ensEMBL2id, ensembl_gene_id == gene)$external_gene_name
    
    tmp.TC.sum          <- ddply(tmp.TC, c("Condition", "Time"), summarise, mean = median(count))
    tmp.TC.sum$gene     <- gene
    tmp.TC.sum$genename <- subset(ensEMBL2id, ensembl_gene_id == gene)$external_gene_name
    
    BigTable         <- rbind(BigTable, tmp.TC)
    BigTable.summary <- rbind(BigTable.summary, tmp.TC.sum)
    count = count+1 
  }
  BigTable %>%
    mutate(genename = factor(genename, levels = as.factor(GENENAME)))
  BigTable.summary %>%
    mutate(genename = factor(genename, levels = as.factor(GENENAME)))
  BigTable$genename <- factor(BigTable$genename, levels = GENENAME)
  plt.TC <-ggplot(BigTable, aes(x = Time, y = count, color = Condition, group = Condition)) + 
    geom_path(data=BigTable.summary,  aes(x=as.factor(Time), y=mean, group=Condition, colour=Condition))  +
    geom_point(data=BigTable.summary, aes(x=as.factor(Time), y=mean, group=Condition, colour=Condition), size=3, alpha=0.75) + 
    geom_jitter(size=1.5, alpha=0.5, shape=1, width=0.1) + 
    scale_color_manual(name="Treatment", values=c("WT"="royalblue", "KO"="darkorange"), limits=c("WT", "KO")) +
    scale_y_log10() +
    xlab("Week") +
    ylab(bquote(paste('log'['10']*'(Normalised Read Counts)'))) +
    ggtitle(TITLE) + 
    facet_wrap( ~ genename, ncol=Fig.ncol) +
    guides(size=FALSE) +
    theme(text=element_text(size=textsize,family="sans"), plot.title=element_text(size=titlesize,family="sans"),
          axis.text.x=element_text(size=axissize,family="sans"), axis.text.y=element_text(size=axissize,family="sans"),
          strip.background = element_rect(colour="red", fill="whitesmoke"),
          legend.position="bottom")
  
  pdf(paste0(Base.dir, "/", Project, "_Fig.TC.", TITLE, ".pdf"),width=Fig.Width,height=Fig.Height)
  par(bg=NA)
  print( plt.TC )
  dev.off()
  
  return(plt.TC)
}
Base.dir <- "/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/March_2023"
Main_markers <- c("Fgf21", "Spp1", "Timp1", "Ager", "Arg1", "Asns", "Psat1", "Phgdh",
                  "Atf3", "Atf4", "Ddit3", "Ccl2", "Ppp1r15a", "Trib3", "Tnfrsf1b",
                  "Il1b", "Gdf15", "Lgals3", "Nlrp3", "Ifi204", "Col11a1", 
                  "Col1a1", "Col3a1", "Col5a1", "Ccr2", "Cebpd", "Fgf2", "Lcn2")
Main_ensembl <- ensEMBL2id[ensEMBL2id$external_gene_name%in%Main_markers==T, c("ensembl_gene_id", "external_gene_name")]
Main_ensembl <- Main_ensembl[-which(duplicated(Main_ensembl$external_gene_name)==T),] # Ddit3 had two ensemblID, ENSMUSG00000025408 & ENSMUSG00000116429, we use the first one
order.ens    <- factor(Main_ensembl$external_gene_name, levels = as.factor(Main_markers))
Main_ensembl <- as.data.frame(Main_ensembl)[order(order.ens),]
plt.1 <- makeCustomGeneSetBubbleplots(Base.dir, Main_ensembl$ensembl_gene_id, Main_markers, dds_all, ensEMBL2id, 15, 10, 4, 
                                      "Main_GO_Genes", 10, 10, 12)         

message("+----------------Immune cell types decovolutions-------------------------------------------+")
message("+---- Microarray had details of cell type compare to RNASeq with a few main categories-----+")

Immune_cells_expr_matrix <- read_excel("./Results_Tables_Figures/March_2023/Deconvolution_Analysis/Decon_Micro_RNA_April_2023.xlsx",  sheet=2)
colnames(Immune_cells_expr_matrix)[1] <- "external_gene_name"

Immune_cells_expr_matrix_anno  <- unique(merge(Immune_cells_expr_matrix, ensEMBL2id.GO, by = "external_gene_name"))
head(Immune_cells_expr_matrix_anno)
Immune_cells_expr_matrix_anno  <- Immune_cells_expr_matrix_anno[order(Immune_cells_expr_matrix_anno$entrezgene),]
Immune_cells_expr_matrix_anno  <- Immune_cells_expr_matrix_anno[!duplicated(Immune_cells_expr_matrix_anno$external_gene_name),]
rownames(Immune_cells_expr_matrix_anno) <- Immune_cells_expr_matrix_anno$ensembl_gene_id
Immune_cells_expr_matrix_anno  <- Immune_cells_expr_matrix_anno[,-c(1,27:28)]

NormCounts_Immune                 <- normcounts[normcounts$ensembl_gene_id %in% rownames(Immune_cells_expr_matrix_anno),c(2:26)]
rownames(NormCounts_Immune) <- NormCounts_Immune[,1]
NormCounts_Immune <- NormCounts_Immune[,-1]
Immune_cells_expr_matrix_anno2    <- (Immune_cells_expr_matrix_anno[rownames(Immune_cells_expr_matrix_anno) %in% rownames(NormCounts_Immune),])

head(NormCounts_Immune)
head(Immune_cells_expr_matrix_anno2)

colnames(NormCounts_Immune)

library(DeconRNASeq)
## due to the limitation of the software can only deal with max 10 cell type each time.

Decon_results.1_10    <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), 
                                     signatures=as.data.frame(Immune_cells_expr_matrix_anno2[,c(1:10)]), 
                                     checksig = F, known.prop = FALSE, use.scale = TRUE, fig = TRUE)
Decon_results.1.10_df <- as.data.frame(Decon_results.1_10$out.all)
rownames(Decon_results.1.10_df) <- colnames(NormCounts_Immune)

Decon_results.11_20    <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), 
                                      signatures=as.data.frame(Immune_cells_expr_matrix_anno2[,c(11:20)]), 
                                      checksig = F, known.prop = FALSE, use.scale = TRUE, fig = TRUE)
Decon_results.11.20_df <- as.data.frame(Decon_results.11_20$out.all)
rownames(Decon_results.11.20_df) <- colnames(NormCounts_Immune)


Decon_results.21_25    <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), 
                                      signatures=as.data.frame(Immune_cells_expr_matrix_anno2[,c(21:25)]), 
                                      checksig = F, known.prop = FALSE, use.scale = TRUE, fig = TRUE)
Decon_results.21.25_df <- as.data.frame(Decon_results.21_25$out.all)
rownames(Decon_results.21.25_df) <- colnames(NormCounts_Immune)


Decon_results_df <- merge(Decon_results.1.10_df, Decon_results.11.20_df,by="row.names")
rownames(Decon_results_df) <- Decon_results_df$Row.names
Decon_results_df$Row.names = NULL

Decon_results_df <- merge(Decon_results_df,  Decon_results.21.25_df,by="row.names",all.x=TRUE)
rownames(Decon_results_df) <- Decon_results_df$Row.names
Decon_results_df$Row.names = NULL


Decon_results_df$condition <- rownames(Decon_results_df)
Decon_results_df$condition <- gsub("_REP[1234]", "", Decon_results_df$condition)
head(Decon_results_df)
colnames(Decon_results_df)[23] <- "NK Actived"


message("+-------------------------------------------------------------------------------+")
message("+ DeconRNA-Seq : Stacked bar plot                                               +")
message("+-------------------------------------------------------------------------------+")

CellTypesList <- c('Plasma.Cells', 'Monocyte', 'Mast.Cells', 'Neutrophil.Cells', 'Eosinophil.Cells',   
                   'NK.Resting', 'NK.Actived', 'M0.Macrophage', 'M1.Macrophage', 'M2.Macrophage',
                   'DC.Immature', 'DC.Actived', 'B.Cells.Memory', 'B.Cells.Naive',  
                   'GammaDelta.T.Cells', 'T.Cells.CD8.Actived', 'T.Cells.CD8.Naive', 
                   'T.Cells.CD8.Memory',  'Treg.Cells',  'T.Cells.CD4.Memory', 'T.Cells.CD4.Naive', 
                   'T.Cells.CD4.Follicular', 'Th1.Cells', 'Th17.Cells', 'Th2.Cells' )
CellTypesList <- gsub("\\.", " ", CellTypesList)
CellTypesList2 <- CellTypesList[-c(13, 22)] ## There is no "B.Cells.Memory" and "T.Cells.CD4.Follicular"
CellTypesList3 <- factor(CellTypesList2[c(2, 8, 9, 10, 11:12, 3, 4, 5, 6, 7, 1, 13, 20, 19, 16, 15, 17, 18, 14,21, 23, 22 )],
                         levels = CellTypesList2[c(2, 8, 9, 10, 11, 12,3, 4, 5, 6, 7, 1, 13, 20, 19, 16, 15, 17, 18, 14,21, 23, 22  )])
 
fill   <- c("firebrick4", "firebrick1", "tomato2", "orangered", "orange3", "orange1", "gold1", "yellow1", "greenyellow",
          "chartreuse3",  "cyan2", "cyan4", "lightskyblue",  "deepskyblue", "dodgerblue", "steelblue2", 
          "royalblue1", "blue3", "navyblue","slateblue2", "darkorchid3", "magenta2", "deeppink2")
fillgr <- c("#990000", "#CC3300", "#FF3300", "#FF6600", 
            "#996600", "#CC9966",
            "#339900", "#66FF33", "#CCFF33",
            "#009999","#33FFFF", 
            "#990099", "#CC66FF",
            "#000066", "#000099", "#003399", "#0033FF", "#3333FF", "#3399FF", "#0099FF", 
            "#006699","#66CCFF", "#99CCFF")
            

Decon_results_df2 <- Decon_results_df[, -which(colnames(Decon_results_df)%in%CellTypesList[c(13,22)]==T) ]

Decon_results_df.means      <- aggregate(Decon_results_df2[,-24], list(Decon_results_df2$condition), mean)
Decon_results_df.means.melt <- melt(Decon_results_df.means)

Decon_results_df.vars      <- aggregate(Decon_results_df2[, -24], list(Decon_results_df2$condition), var)
Decon_results_df.vars.melt <- melt(Decon_results_df.vars)

Decon_results_df.means.melt$vars     <- Decon_results_df.vars.melt$value

# Calculate a normalised variance
Decon_results_df.means.melt$normvars <- Decon_results_df.means.melt$vars/max(Decon_results_df.means.melt$vars)

# Tidy up the table a bit for plotting
colnames(Decon_results_df.means.melt) <- c("Group", "Immune.Cell", "mean", "variance", "normvars")

Decon_results_df.means.melt$Immune.Cell <- factor(Decon_results_df.means.melt$Immune.Cell, levels =  CellTypesList3 )
head(Decon_results_df.means.melt)



Decon_results_df.means.melt2 <- Decon_results_df.means.melt
Decon_results_df.means.melt2$condition <- substr(Decon_results_df.means.melt$Group, 1,2)
Decon_results_df.means.melt2$condition <- factor(Decon_results_df.means.melt2$condition, levels=c("WT", "KO"))
Decon_results_df.means.melt2$time <- substr(Decon_results_df.means.melt$Group, 5, 5)
Decon_results_df.means.melt2$time      <- factor(Decon_results_df.means.melt2$time, levels=c("0", "2", "4", "7"))
Decon.stacked.bar <- ggplot(data = Decon_results_df.means.melt2, aes(y = mean, x = time, fill = Immune.Cell)) + 
  geom_bar( colour = "black", stat="identity", position = "stack") +  
  scale_fill_manual(values=fillgr) +
  # scale_x_discrete(name="", limits=c("WT_T0w", "WT_T2w", "WT_T4w", "WT_T7w","KO_T0w", "KO_T2w","KO_T4w", "KO_T7w"),
  #                 labels=c("0", "2","4","7", "0", "2","4","7")) +
  scale_y_continuous(breaks=c(0,1,2,3), labels=c(0,0.25,0.75,1))+
  theme(legend.title = element_blank(), legend.key.height = unit(3.0, 'cm')) + 
  labs(x="", y="Fraction of cells") + 
  facet_grid(~condition) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=14, face="bold"),
        axis.title.y=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face= "bold"))+
  theme(strip.text.x = element_text(size = 14, face = "bold")) +
  guides(fill = guide_legend(ncol = 1)) 
# revised with Jane on 24/04/2023 
pdf(paste("./Results_Tables_Figures/Paper_v1_April_2023/Figxx-DESeq2_ImmuneCellTypeDeconvolution_stacked_barplot_May_2023.pdf", sep=""), width=10,height=8, onefile=FALSE)
par(bg=NA)
Decon.stacked.bar
dev.off()  

#######################################################
pdf(paste("./Results_Tables_Figures/March_2023/Deconvolution_Analysis/", Project, "-Figxx_DESeq2_ImmuneCellTypeDeconvolution_stacked_barplot.pdf", sep=""), width=10,height=8, onefile=FALSE)
par(bg=NA)
Decon.stacked.bar
dev.off()
#######################################################

message("+-------------------------------------------------------------------------------+")
message("+ GeneOntology Pathway Heatmap with Genes                                       +")
message("+-------------------------------------------------------------------------------+")

pathwayNames <- c("ATF4 activates \nin response to \nER stress & ISR",
                  "Fibrosis",
                  "Pyruvate metabolism \nand TCA cycle", 
                  "Serine biosynthesis \nand one carbon cycle", 
                  "p53 pathway",
                  "Intrinsic apoptotic \nsignaling pathway \nin response to \nDNA damage",
                  "Inflammation",
                  "Cytokine-cytokine \nreceptor interaction",
                  "Response to \noxidative stress")
pathGenesList <- list(c("Trib3", "Atf3", "Asns", "Ddit3", "Ppp1r15a", "Chac1", "Sesn2", "Atf4", "Ccl2"),
                      c("Col11a1", "Col16a1", "Lgals3", "Ccn2", "Tnfrsf11b", "Mmp24", "Lox", "Fn1"),
                      c("Sdh", "Fh1", "Mpc1", "Pdp1", "Pdk1"),
                      c("Arg1", "Mthfd2", "Aldh1l2", "Phgdh", "Psat1", "Psph"),
                      c("Serpine1", "Fas", "Pidd1", "Gadd45a", "Cdkn1a", "Trp53inp1", "Bax"),
                      c("Hmox1", "Ddit4", "Bcl2a1a", "Muc1", "Tnfrsf1b", "Sfn"),
                      c("Nlrp3", "Il18", "Nfkb2", "Phlda3", "Pycard", "Casp1", "Ifi204", "Il1b"), 
                      c("Cnn2", "Gdf15", "Ccl8", "Cxcl10", "Ccr2", "Tnfsf18", "Spp1"),
                      c("Ide", "Tnfaip3", "Acox2", "Trem2", "Fos", "Lcn2", "Cebpb"))
							
degs_wk4 <- read.csv(RES.files1[3], header = T)[,c("external_gene_name", "log2FoldChange", "padj")]
degs_wk7 <- read.csv(RES.files1[4], header = T)[,c("external_gene_name", "log2FoldChange", "padj")]
degs_wk20 <- read.csv(RESWT.files[1], header = T)[,c("external_gene_name", "log2FoldChange", "padj")] 
degs_wk40 <- read.csv(RESWT.files[2], header = T)[,c("external_gene_name", "log2FoldChange", "padj")]
degs_wk70 <- read.csv(RESWT.files[3], header = T)[,c("external_gene_name", "log2FoldChange", "padj")] 

pathwk4ko <- lapply(pathGenesList, function(x) degs_wk4[degs_wk4$external_gene_name%in%x,])
pathwk7ko <- lapply(pathGenesList, function(x) degs_wk7[degs_wk7$external_gene_name%in%x,])

pathwk2wt <- lapply(pathGenesList, function(x) degs_wk20[degs_wk20$external_gene_name%in%x,])
pathwk4wt <- lapply(pathGenesList, function(x) degs_wk40[degs_wk40$external_gene_name%in%x,])
pathwk7wt <- lapply(pathGenesList, function(x) degs_wk70[degs_wk70$external_gene_name%in%x,])


pldList <- list()
for(i in 1:length(pathGenesList)){
  pld <- Reduce(function(df1, df2) merge(df1, df2, by = "external_gene_name", all=T), 
                         list(pathwk4ko[[i]], pathwk7ko[[i]], pathwk2wt[[i]], 
                              pathwk4wt[[i]], pathwk7wt[[i]]))
  colnames(pld) <- c("GeneName", "w4KO_l2fc", "w4KO_padj", "w7KO_l2fc", "w7KO_padj",
                      "w2WT_l2fc", "w2WT_padj","w4WT_l2fc", "w4WT_padj", "w7WT_l2fc", "w7WT_padj")
  rownames(pld) <- pld$GeneName
  pld.l2fc <- pld[,c(1, 2, 4, 6, 8, 10)]; pld.padj <- pld[,c(1, 3, 5, 7, 9, 11)]
  pld.l2fcm <- melt(pld.l2fc, "GeneName"); pld.padjm <- melt(pld.padj, "GeneName"); 
  colnames(pld.l2fcm) <- c("GeneName", "Group", "L2FC")
  colnames(pld.padjm) <- c("GeneName", "Group", "padj")
  pldList[[i]] <- cbind(pld.l2fcm, pld.padjm[,3])
  colnames(pldList[[i]]) <- c("GeneName", "Group","L2FC", "padj")
  pldList[[i]][,2] <- gsub("_l2fc", "", pldList[[i]][,2])
  pldList
}

pldListHeat <- list()
for(i in 1:length(pathGenesList)){
  pld <- Reduce(function(df1, df2) merge(df1, df2, by = "external_gene_name", all=T), 
                list(pathwk4ko[[i]], pathwk7ko[[i]], pathwk2wt[[i]], 
                     pathwk4wt[[i]], pathwk7wt[[i]]))
  colnames(pld) <- c("GeneName", "w4KO_l2fc", "w4KO_padj", "w7KO_l2fc", "w7KO_padj",
                     "w2WT_l2fc", "w2WT_padj","w4WT_l2fc", "w4WT_padj", "w7WT_l2fc", "w7WT_padj")
  rownames(pld) <- pld$GeneName
  
  pldListHeat[[i]] <- pld[,c(1,2,4,6,8,10)]
  pldListHeat[[i]]$Term <- pathwayNames[i]
  pldListHeat
}


colsf1 = colorRamp2( c(-2,-0.6,0,0.6,5), c("darkgreen", "green3", "gray95", "slateblue", "darkorchid4"), space = "RGB") 

pltall <- data.table::rbindlist(pldListHeat)
pltall <- as.data.frame(pltall)
rownames(pltall) <- pltall$GeneName
pltallHeat <- ComplexHeatmap:: Heatmap(as.matrix(pltall[,2:6]), na_col ="black",
                           col = colsf1, name = "pathheat",  
                           row_title = "", column_title = NULL, 
                           show_row_names = TRUE, show_column_names = TRUE,
                           column_labels = c(4, 7, 2, 4, 7), 
                           heatmap_legend_param = list(title = "log2FoldChange", 
                                                       legend_height = unit(6, "cm"), 
                                                       title_position = "leftcenter-rot"), 
                           cluster_columns = FALSE, row_names_side ="right", 
                           row_split = pltall$Term, cluster_rows = TRUE, 
                           column_split = c("KOvsWT", "KOvsWT", "WT2/4/7vsWT0","WT2/4/7vsWT0","WT2/4/7vsWT0"),
                           column_names_rot = 0, column_names_side = "bottom", 
                           column_title_rot = 0, column_names_centered = TRUE,
                           column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                           row_names_gp = gpar(fontsize = 10, fontface = "bold.italic"),
                           column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                           row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                           row_title_rot = 0, cluster_row_slices = FALSE, 
                           show_parent_dend_line = FALSE,
                           width = ncol(pltall)*unit(5, "mm"), 
                           height = nrow(pltall)*unit(4, "mm"),
                           row_gap = unit(2, "mm"), column_gap = unit(3, "mm"))
  
pdf(paste("./Results_Tables_Figures/Paper_v1_April_2023/Figxx-SelPathways_genes_heatmap_May_2023.pdf", sep=""), width=10,height=14, onefile=FALSE)
par(bg=NA)
print(pltallHeat)
dev.off()  

message("+--- Review of the BHJ regarding to the Figure 2c and SFig2g ----------------+")

## The reviewer asked whether these pathways are up/dw regulated
## Solutions: due to our heatmap is merging some pathways relating to the interests, and it is not quite straightforward for answering.
## I will gather these genes in heatmap and check the overlap list with pathways and get the up/dw question based on the genes up/dw direction.
library(rlist)
goFile <- "/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/Paper_v1_April_2023/Supplementary_Data_May_2023.xlsx"

koresfiles <- paste0("/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/March_2023/CAD_jcg23_0001-Resall_lfcShrink_",
                     c("4w", "7w"), "_summaryTable_March_2023.csv")
pathGenesList <- list(c("Trib3", "Atf3", "Asns", "Ddit3", "Ppp1r15a", "Chac1", "Sesn2", "Atf4", "Ccl2"),
                      c("Col11a1", "Col16a1", "Lgals3", "Ccn2", "Tnfrsf11b", "Mmp24", "Lox", "Fn1"),
                      c("Sdh", "Fh1", "Mpc1", "Pdp1", "Pdk1"),
                      c("Arg1", "Mthfd2", "Aldh1l2", "Phgdh", "Psat1", "Psph"),
                      c("Serpine1", "Fas", "Pidd1", "Gadd45a", "Cdkn1a", "Trp53inp1", "Bax"),
                      c("Hmox1", "Ddit4", "Bcl2a1a", "Muc1", "Tnfrsf1b", "Sfn"),
                      c("Nlrp3", "Il18", "Nfkb2", "Phlda3", "Pycard", "Casp1", "Ifi204", "Il1b"), 
                      c("Cnn2", "Gdf15", "Ccl8", "Cxcl10", "Ccr2", "Tnfsf18", "Spp1"),
                      c("Ide", "Tnfaip3", "Acox2", "Trem2", "Fos", "Lcn2", "Cebpb"))
pathNamesList <- c("ER/ISR", "Cytokine", "EMO", "Inflam", "DNAdamage", "P53", "TCA", "Oxidative", "Serine/OneCarbon")
gthrKO_fn <- function(inputfile, resfiles=koresfiles,  sindex,
                      konames, pgenes, pgnames, outfile){
  dat  <- openxlsx::read.xlsx(inputfile, sheet = sindex, startRow = 2)
  res1 <- read.csv(resfiles[1], header = T)
  res2 <- read.csv(resfiles[2], header = T)
  GeneLists1 <- lapply(dat[dat$Cluster==konames[1],"geneID"], function(x) strsplit(x, split = "/"))
  GeneLists2 <- lapply(dat[dat$Cluster==konames[2],"geneID"], function(x) strsplit(x, split = "/"))
  GeneLists <- c(GeneLists1, GeneLists2)
  
  geover1 <- lapply(GeneLists1, function(x) res1[res1$external_gene_name%in%unlist(x), "log2FoldChange"])
  geover1_sign <-  unlist(lapply(geover1, function(x) sign(sum(x>0) - sum(x<0))))
  
  geover2 <- lapply(GeneLists2, function(x) res2[res2$external_gene_name%in%unlist(x), "log2FoldChange"])
  geover2_sign <- unlist(lapply(geover2, function(x) sign(sum(x>0) - sum(x<0))))
  
  geover_sign <- c(geover1_sign, geover2_sign)
  geover_sign <- ifelse(geover_sign == "-1", "Down", geover_sign)
  geover_sign <- ifelse(geover_sign == "1", "Up", geover_sign)
  geover_sign <- ifelse(geover_sign == "0", "Equal", geover_sign)
  
  dat$PathwayDirection <- geover_sign
  pgoverall <- list()
  for(i in 1:length(pgenes)){
    pgoverall[[i]] <- unlist(lapply(GeneLists, function(x) paste0(unlist(x)[unlist(x)%in%pgenes[[i]]], collapse="/")))
  }
  pgoverdat <- list.cbind(pgoverall)
  colnames(pgoverdat) <- pgnames
  ndat <- cbind(dat, pgoverdat)
  openxlsx::write.xlsx(ndat, file = outfile)
  gc()
}

gowk47KO <- gthrKO_fn(inputfile=goFile, resfiles=koresfiles,  sindex = 10,
                      konames=c("wk4", "wk7"), pgenes=pathGenesList, pgnames=pathNamesList, 
                      outfile="/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/March_2023/Review_EHJ_GO_KO_addDirection_10102023.xlsx")

keggwk47KO <- gthrKO_fn(inputfile=goFile, resfiles=koresfiles,  sindex = 12,
                         konames=c("wk4", "wk7"), pgenes=pathGenesList, pgnames=pathNamesList, 
                         outfile="/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/March_2023/Review_EHJ_KEGG_KO_addDirection_10102023.xlsx")
gthrWT_fn <- function(inputfile, resfiles=wtresfiles,  sindex,
                      wtnames, pgenes, pgnames, outfile){
  dat  <- openxlsx::read.xlsx(inputfile, sheet = sindex, startRow = 2)
  res1 <- read.csv(resfiles[1], header = T)
  res2 <- read.csv(resfiles[2], header = T)
  res3 <- read.csv(resfiles[3], header = T)
  GeneLists1 <- lapply(dat[dat$Cluster==wtnames[1],"geneID"], function(x) strsplit(x, split = "/"))
  GeneLists2 <- lapply(dat[dat$Cluster==wtnames[2],"geneID"], function(x) strsplit(x, split = "/"))
  GeneLists3 <- lapply(dat[dat$Cluster==wtnames[3],"geneID"], function(x) strsplit(x, split = "/"))
  GeneLists <- c(GeneLists1, GeneLists2, GeneLists3)
  
  geover1 <- lapply(GeneLists1, function(x) res1[res1$external_gene_name%in%unlist(x), "log2FoldChange"])
  geover1_sign <-  unlist(lapply(geover1, function(x) sign(sum(x>0) - sum(x<0))))
  
  geover2 <- lapply(GeneLists2, function(x) res2[res2$external_gene_name%in%unlist(x), "log2FoldChange"])
  geover2_sign <- unlist(lapply(geover2, function(x) sign(sum(x>0) - sum(x<0))))
  
  geover3 <- lapply(GeneLists3, function(x) res3[res3$external_gene_name%in%unlist(x), "log2FoldChange"])
  geover3_sign <- unlist(lapply(geover3, function(x) sign(sum(x>0) - sum(x<0))))
  
  geover_sign <- c(geover1_sign, geover2_sign, geover3_sign)
  geover_sign <- ifelse(geover_sign == "-1", "Down", geover_sign)
  geover_sign <- ifelse(geover_sign == "1", "Up", geover_sign)
  geover_sign <- ifelse(geover_sign == "0", "Equal", geover_sign)
  
  dat$PathwayDirection <- geover_sign
  pgoverall <- list()
  for(i in 1:length(pgenes)){
    pgoverall[[i]] <- unlist(lapply(GeneLists, function(x) paste0(unlist(x)[unlist(x)%in%pgenes[[i]]], collapse="/")))
  }
  pgoverdat <- list.cbind(pgoverall)
  colnames(pgoverdat) <- pgnames
  ndat <- cbind(dat, pgoverdat)
  openxlsx::write.xlsx(ndat, file = outfile)
  gc()
}
wtresfiles <- paste0("/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/March_2023/CAD_jcg23_0001-Resall_WT_lfcShrink_",
                     c("2wvs0w", "4wvs0w", "7wvs0w"), "_summaryTable_March_2023.csv")
gowk247WT <- gthrWT_fn(inputfile=goFile,  resfiles=wtresfiles,  sindex = 13,
                      wtnames=c("wk2", "wk4", "wk7"), pgenes=pathGenesList, pgnames=pathNamesList, 
                      outfile="/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/March_2023/Review_EHJ_GO_WT_addDirection_10102023.xlsx")

keggwk247WT <- gthrWT_fn(inputfile=goFile, resfiles=wtresfiles,  sindex = 15,
                        wtnames=c("wk2", "wk4", "wk7"), pgenes=pathGenesList, pgnames=pathNamesList, 
                        outfile="/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/March_2023/Review_EHJ_KEGG_WT_addDirection_10102023.xlsx")


## Reactome is used orthology human genes, need to call the data .

load(file = "./Results_Tables_Figures/March_2023/Gene_Ontology/Reactome_inputortho_KO_wk4_wk7_data.RData")
kowk4 <- wk4_geneList; kowk7 <- wk7_geneList
load(file = "./Results_Tables_Figures/March_2023/Gene_Ontology/WT_wk247vswk0_reactome_input_bkg_geneList_March_2023.RData")
wtwk2 <- WT.orthos[[1]]
wtwk4 <- WT.orthos[[2]]
wtwk7 <- WT.orthos[[3]]

koreact <- openxlsx::read.xlsx(goFile, sheet = 11, startRow = 2)
wtreact <- openxlsx::read.xlsx(goFile, sheet = 14, startRow = 2)

reacthrKO_fn <- function(inputfile, resfile1=kowk4, resfile2=kowk7,  sindex,
                      konames, pgenes, pgnames, outfile){
  dat  <- openxlsx::read.xlsx(inputfile, sheet = sindex, startRow = 2)
  
  GeneLists1 <- lapply(dat[dat$Cluster==konames[1],"geneID"], function(x) strsplit(x, split = "/"))
  GeneLists2 <- lapply(dat[dat$Cluster==konames[2],"geneID"], function(x) strsplit(x, split = "/"))

  geover1 <- lapply(GeneLists1, function(x) unique(resfile1[resfile1$Human.gene.name%in%unlist(x), c("SYMBOL", "L2FC")]))
  geover1_sign <-  unlist(lapply(geover1, function(x) sign(sum(x$L2FC>0) - sum(x$L<0))))
  geneList1 <- lapply(geover1, function(x) x$SYMBOL)
  geover2 <- lapply(GeneLists2, function(x) unique(resfile2[resfile2$Human.gene.name%in%unlist(x), c("SYMBOL", "L2FC")]))
  geover2_sign <- unlist(lapply(geover2, function(x) sign(sum(x$L2FC>0) - sum(x$L2FC<0))))
  geneList2 <- lapply(geover2, function(x) x$SYMBOL)
  GeneLists <- c(geneList1, geneList2)
  
  geover_sign <- c(geover1_sign, geover2_sign)
  geover_sign <- ifelse(geover_sign == "-1", "Down", geover_sign)
  geover_sign <- ifelse(geover_sign == "1", "Up", geover_sign)
  geover_sign <- ifelse(geover_sign == "0", "Equal", geover_sign)
  
  dat$PathwayDirection <- geover_sign
  pgoverall <- list()
  for(i in 1:length(pgenes)){
    pgoverall[[i]] <- unlist(lapply(GeneLists, function(x) paste0(unlist(x)[unlist(x)%in%pgenes[[i]]], collapse="/")))
  }
  pgoverdat <- list.cbind(pgoverall)
  colnames(pgoverdat) <- pgnames
  ndat <- cbind(dat, pgoverdat)
  openxlsx::write.xlsx(ndat, file = outfile)
  gc()
}
reacthrKO_fn(inputfile=goFile, resfile1=kowk4, resfile2=kowk7,  sindex=11, konames=c("wk4", "wk7"), pgenes=pathGenesList, 
             pgnames=pathNamesList, 
             outfile="/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/March_2023/Review_EHJ_Reactome_KO_addDirection_10102023.xlsx")

reacthrWT_fn <- function(inputfile, pathG, resfile1=wtwk2, resfile2=wtwk4, resfile3=wtwk7, sindex,
                         wtnames, pgenes, pgnames, outfile){
  dat  <- openxlsx::read.xlsx(inputfile, sheet = sindex, startRow = 2)
  
  GeneLists1 <- lapply(dat[dat$Cluster==wtnames[1],"geneID"], function(x) strsplit(x, split = "/"))
  GeneLists2 <- lapply(dat[dat$Cluster==wtnames[2],"geneID"], function(x) strsplit(x, split = "/"))
  GeneLists3 <- lapply(dat[dat$Cluster==wtnames[3],"geneID"], function(x) strsplit(x, split = "/"))

  geover1 <- lapply(GeneLists1, function(x) unique(resfile1[resfile1$Human.gene.name%in%unlist(x), c("SYMBOL", "L2FC")]))
  geover1_sign <-  unlist(lapply(geover1, function(x) sign(sum(x$L2FC>0) - sum(x$L<0))))
  geneList1 <- lapply(geover1, function(x) x$SYMBOL)
  geover2 <- lapply(GeneLists2, function(x) unique(resfile2[resfile2$Human.gene.name%in%unlist(x), c("SYMBOL", "L2FC")]))
  geover2_sign <- unlist(lapply(geover2, function(x) sign(sum(x$L2FC>0) - sum(x$L2FC<0))))
  geneList2 <- lapply(geover2, function(x) x$SYMBOL)
  geover3 <- lapply(GeneLists3, function(x) unique(resfile3[resfile3$Human.gene.name%in%unlist(x), c("SYMBOL", "L2FC")]))
  geover3_sign <- unlist(lapply(geover3, function(x) sign(sum(x$L2FC>0) - sum(x$L2FC<0))))
  geneList3 <- lapply(geover3, function(x) x$SYMBOL)
  GeneLists <- c(geneList1, geneList2, geneList3)
  
  geover_sign <- c(geover1_sign, geover2_sign, geover3_sign)
  geover_sign <- ifelse(geover_sign == "-1", "Down", geover_sign)
  geover_sign <- ifelse(geover_sign == "1", "Up", geover_sign)
  geover_sign <- ifelse(geover_sign == "0", "Equal", geover_sign)
  
  dat$PathwayDirection <- geover_sign
  pgoverall <- list()
  for(i in 1:length(pgenes)){
    pgoverall[[i]] <- unlist(lapply(GeneLists, function(x) paste0(unlist(x)[unlist(x)%in%pgenes[[i]]], collapse="/")))
  }
  pgoverdat <- list.cbind(pgoverall)
  colnames(pgoverdat) <- pgnames
  ndat <- cbind(dat, pgoverdat)
  openxlsx::write.xlsx(ndat, file = outfile)
  gc()
}
reacthrWT_fn(inputfile=goFile, resfile1=wtwk2, resfile2=wtwk4,  resfile3=wtwk7, sindex=14, wtnames=c("wk2", "wk4", "wk7"), 
             pgenes=pathGenesList, pgnames=pathNamesList, 
             outfile="/Users/xz289/Documents/Jane_Goodall/CAD_jcg23_0001/Results_Tables_Figures/March_2023/Review_EHJ_Reactome_WT_addDirection_10102023.xlsx")
##-------------------------------Finish 19/10/2023 ----------------------------------------------##
