library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(gplots)
library(ggpubr)
library(tidyverse)
library(DOSE)
library(KEGG.db)
library(pathview)
library(graphics)

# creating dataframes

df1 = read.csv("sig_Low_AL_DESeq2.csv")
head(df1)

df2 = read.csv("sig_High_AL_DESeq2.csv")
head(df2)

df3 = read.csv("sig_Low_S1133_DESeq2.csv")
head(df3)

df4 = read.csv("sig_High_S1133_DESeq2.csv")
head(df4)

# Log2FoldChange

df1_gene_list <- df1$log2FoldChange
head(df1_gene_list)

df2_gene_list <- df2$log2FoldChange
head(df2_gene_list)

df3_gene_list <- df3$log2FoldChange
head(df3_gene_list)

df4_gene_list <- df4$log2FoldChange
head(df4_gene_list)

# name the vector
names(df1_gene_list) <- df1$X
names(df2_gene_list) <- df2$X
names(df3_gene_list) <- df3$X
names(df4_gene_list) <- df4$X


# omit any NA values 
df1_gene_list<-na.omit(df1_gene_list)
df2_gene_list<-na.omit(df2_gene_list)
df3_gene_list<-na.omit(df3_gene_list)
df4_gene_list<-na.omit(df4_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
df1_gene_list = sort(df1_gene_list, decreasing = TRUE)
head(df1_gene_list)

df2_gene_list = sort(df2_gene_list, decreasing = TRUE)
head(df2_gene_list)

df3_gene_list = sort(df3_gene_list, decreasing = TRUE)
head(df3_gene_list)

df4_gene_list = sort(df4_gene_list, decreasing = TRUE)
head(df4_gene_list)



# compare groups
upregulated_genes_Low_AL <- subset(df1, log2FoldChange > 1)
downregulated_genes_Low_AL <- subset(df1, log2FoldChange < -1)
upregulated_genes_High_AL <- subset(df2, log2FoldChange > 1)
downregulated_genes_High_AL <- subset(df2, log2FoldChange < -1)
upregulated_genes_Low_S1133 <- subset(df3, log2FoldChange > 1)
downregulated_genes_Low_S1133 <- subset(df3, log2FoldChange < -1)
upregulated_genes_High_S1133 <- subset(df4, log2FoldChange > 1)
downregulated_genes_High_S1133 <- subset(df4, log2FoldChange < -1)

names(upregulated_genes_High_AL)

ck <- compareCluster(list(
  Low_AL_up=upregulated_genes_Low_AL$X, Low_AL_down=downregulated_genes_Low_AL$X,
  High_AL_up=upregulated_genes_High_AL$X, High_AL_down=downregulated_genes_High_AL$X,
  Low_S1133_up=upregulated_genes_Low_S1133$X, Low_S1133_down=downregulated_genes_Low_S1133$X,
  High_S1133_up=upregulated_genes_High_S1133$X, High_S1133_down=downregulated_genes_High_S1133$X),
  fun = "enrichKEGG", organism = "ko",
  pvalueCutoff = 0.05)

#High_AL_up=names(upregulated_genes_High_AL), High_AL_down=names(downregulated_genes_High_AL),
#Low_S1133_up=names(upregulated_genes_Low_S1133), Low_S1133_down=names(downregulated_genes_Low_S1133),
#High_S1133_up=names(upregulated_genes_High_S1133), High_S1133_down=names(downregulated_genes_High_S1133)),

dotplot(ck)