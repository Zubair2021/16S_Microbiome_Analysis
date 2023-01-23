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

# we want the log2 fold change 

df1_gene_list <- df1$log2FoldChange
head(df1_gene_list)

df2_gene_list <- df2$log2FoldChange
head(df2_gene_list)

df3_gene_list <- df3$log2FoldChange
head(df3_gene_list)

df4_gene_list <- df4$log2FoldChange
head(df4_gene_list)


# Define upregulated



# name the vector
names(df1_gene_list) <- df1$X
names(df2_gene_list) <- df2$X
names(df3_gene_list) <- df3$X
names(df4_gene_list) <- df4$X

# KO_names <- names(original_gene_list)
# head(KO_names)

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

# gene names
#gene_names <- names(gene_list)
#head(gene_names)

# compare groups

ck <- compareCluster(list(Low_AL=names(df1_gene_list), High_AL=names(df2_gene_list),
                    Low_S1133=names(df3_gene_list), High_S1133=names(df4_gene_list)),
                    fun = "enrichKEGG", organism = "ko",
                    pvalueCutoff = 0.05)

dotplot(ck, showCategory = 5)
cnetplot(ck)
ridgeplot(ck)

#compare modules in groups

Mck <- compareCluster(list(Low_AL=names(df1_gene_list), High_AL=names(df2_gene_list),
                          Low_S1133=names(df3_gene_list), High_S1133=names(df4_gene_list)),
                     fun = "enrichMKEGG", organism = "ko",
                     pvalueCutoff = 0.05)
dotplot(Mck, showCategory = 20)



cnetplot(Mck)
