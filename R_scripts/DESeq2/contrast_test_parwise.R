library("phyloseq"); packageVersion("phyloseq")
library(tidyverse)
library(gplots)
library("DESeq2")
packageVersion("DESeq2")
library("dplyr")


# Importing data
data <- read.table("sorted_abun.tsv", header = T)
head(data)
ncol(data)

# Rounding off
data <- data %>% mutate_if(is.numeric, ~round(., 0))
head(data)

# metadata
metadata <- read.csv("metadata-R.tsv", sep = "\t", na.strings = c("NA", ""), header = T)
metadata <- na.omit(metadata)
metadata
nrow(metadata) # make sure the ncol of data is the same as nrow of metadata



# create the DESeq dataset

dds <- DESeqDataSetFromMatrix(countData = data, 
                                              colData = metadata,
                                              design = ~ Treatment)

# test contrast 

#Deseq_contrast_data$treatment <- relevel(Deseq_contrast_data$treatment, ref = "control")

Deseq_con_MLE <- DESeq(Deseq_contrast_data, modelMatrixType="expanded", betaPrior=T)
Deseq_con_MLE

# Checking group names
resultsNames(Deseq_con_MLE)

# defining contrasts

contrast1 = c("Treatment", "Low_AL", "Control")
contrast2 = c("Treatment", "High_AL", "Control")
contrast3 = c("Treatment", "Low_S1133", "Control")
contrast4 = c("Treatment", "High_S1133", "Control")

# running contrasts
contrast1_res_MLE <- results(Deseq_con_MLE, contrast=contrast1, alpha = 0.05)
contrast2_res_MLE <- results(Deseq_con_MLE, contrast=contrast2, alpha = 0.05)
contrast3_res_MLE <- results(Deseq_con_MLE, contrast=contrast3, alpha = 0.05)
contrast4_res_MLE <- results(Deseq_con_MLE, contrast=contrast4, alpha = 0.05)


# writing csv of all results
write.csv(as.matrix(contrast1_res_MLE), file = "all_Low_AL_DESeq2.csv", row.names = T)
write.csv(as.matrix(contrast2_res_MLE), file = "all_High_AL_DESeq2.csv", row.names = T)
write.csv(as.matrix(contrast3_res_MLE), file = "all_Low_S1133_DESeq2.csv", row.names = T)
write.csv(as.matrix(contrast4_res_MLE), file = "all_High_S1133_DESeq2.csv", row.names = T)

# significant results
as.data.frame(contrast1_res_MLE) %>% 
  filter(padj<0.05) -> sig_contrast1_MLE

as.data.frame(contrast2_res_MLE) %>% 
  filter(padj<0.05) -> sig_contrast2_MLE

as.data.frame(contrast3_res_MLE) %>% 
  filter(padj<0.05) -> sig_contrast3_MLE

as.data.frame(contrast4_res_MLE) %>% 
  filter(padj<0.05) -> sig_contrast4_MLE

# writing csv of significant results from contrasts
write.csv(as.matrix(sig_contrast1_MLE), file = "sig_Low_AL_DESeq2.csv", row.names = T)
write.csv(as.matrix(sig_contrast2_MLE), file = "sig_High_AL_DESeq2.csv", row.names = T)
write.csv(as.matrix(sig_contrast3_MLE), file = "sig_Low_S1133_DESeq2.csv", row.names = T)
write.csv(as.matrix(sig_contrast4_MLE), file = "sig_High_S1133_DESeq2.csv", row.names = T)

