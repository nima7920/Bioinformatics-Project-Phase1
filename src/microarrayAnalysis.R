setwd("E:\\University\\Bioinformatics\\Project\\Bioinformatics-Project-Phase1")
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(plyr)
library(Rtsne)

# loading dataset
gset <- getGEO("GSE48558", GSEMatrix =TRUE, getGPL=FALSE,destdir = "data/")

gset <- gset[[1]]
matrix <- exprs(gset) 
# extracts the matrix from dataset

# removing useless samples
patients<-matrix[,1:13]
normal_indices=c(41,45,69,71,75,77:80,82,85,86,92,94,96,98,100,104,108,138:144,147,151:170)
normals <-matrix[,normal_indices]
matrix<-cbind(patients,normals)
labels <-c(rep("AML",13),rep("Normal",47))

# normalizaing features
matrix <-normalizeQuantiles(matrix)

# boxplot of each feature
pdf("results/boxplot.pdf",width=170)
boxplot(matrix)
dev.off()

# heatmap of correlation between features
pdf("results/Correlation heatmap.pdf",width=10,height=10)
pheatmap(cor(matrix),labels_row = labels,labels_col = labels)
dev.off()

#### PCA
# pca for genes
pcs <- prcomp(matrix)
x_transformed <- pcs$x
pdf("results/pcs.pdf")
plot(pcs)
plot(x_transformed[,1:2])
dev.off()

# scaling to extract the difference between gene expressions
matrix.scale <-t(scale(t(matrix),scale=F))
pcs <-prcomp(matrix.scale)
x_transformed <- pcs$x
pdf("results/pcs_scaled.pdf")
plot(pcs)
plot(x_transformed[,1:2])
dev.off()

# PCA for samples
# labels should be added , 1:15
pcr <-data.frame(pcs$rotation[,1:3],Group=labels)
pdf("results/samples pca.pdf")
ggplot(pcr,aes(PC1,PC2,color=Group)) + geom_point(size=2) + theme_bw()
dev.off()

#### MDS for samples
distance_matrix <- dist(t(matrix))
mds <-cmdscale(distance_matrix,k=2)
mds <-data.frame(mds[,1:2],Group=labels)
pdf("results/samples_MDS.pdf")
# plot(mds[,1],mds[,2],xlab="MDS1",ylab="MDS2",main="MDS")
ggplot(mds,aes(mds[,1],mds[,2],color=Group))+ geom_point(size=2) + theme_bw()
dev.off()

#### t-SNE for samples
tsne<-Rtsne(t(matrix),perplexity = 15)
tsne<-data.frame(tsne$Y,Group=labels)
pdf("results/samples_tsne.pdf")
ggplot(tsne,aes(tsne[,1],tsne[,2],color=Group))+ geom_point(size=2) + theme_bw()
#plot(tsne$Y)
dev.off()

