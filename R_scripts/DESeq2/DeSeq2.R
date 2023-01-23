library("phyloseq"); packageVersion("phyloseq")
library(tidyverse)
library(gplots)
library("DESeq2")
packageVersion("DESeq2")


# Importing data
data <- read.table("sorted_abun_high_AL.txt", header = T)
head(data)
ncol(data)

metadata <- read.csv("metadata.tsv", sep = "\t", na.strings = c("NA", ""), header = T)
metadata <- na.omit(metadata)
metadata
nrow(metadata)

#metadata <- read.table("metadata-R.tsv", header = T)
#head(metadata)

# Rounding off only numeric columns to non-decimals using dplyr
data <- data %>% mutate_if(is.numeric, ~round(., 0))
head(data)

# ?DESeqDataSet()

# Deseq_data <- DESeqDataSetFromMatrix(countData = data, 
Deseq_data <- DESeqDataSetFromMatrix(countData = data, 
                                     colData = metadata,
                                     design = ~ specific)
summary(Deseq_data)


# Including the sampling_day from abundance and metadata matrices 

Deseq_data2 <- DESeqDataSetFromMatrix(countData = data, 
                                      colData = metadata,
                                      design = ~ specific + sampling_day)
Deseq_data2
summary(Deseq_data2)

# with interactions bw response variables
Deseq_data3 <- DESeqDataSetFromMatrix(countData = data, 
                                      colData = metadata,
                                     design = ~ specific + sampling_day + specific:sampling_day)
Deseq_data3


# reordering factors and specifying the uninfected control 

Deseq_data$specific <- relevel(Deseq_data$specific, ref = "control")
Deseq_data$specific

# data2
Deseq_data2$specific <- relevel(Deseq_data2$specific, ref = "control")
Deseq_data2$specific

# data 3
Deseq_data3$specific <- relevel(Deseq_data3$specific, ref = "control")
Deseq_data3$specific


# MLE

#?DESeq

Deseq_MLE <- DESeq(Deseq_data, modelMatrixType="standard", betaPrior=FALSE)
Deseq_MLE

Deseq2_MLE <- DESeq(Deseq_data2, modelMatrixType="standard", betaPrior=FALSE)
Deseq2_MLE

Deseq3_MLE <- DESeq(Deseq_data3, modelMatrixType="standard", betaPrior=FALSE)
Deseq3_MLE

# MLE Results
a = 0.05

#?results
res_MLE <- results(Deseq_MLE, name="specific_HighAL_vs_control", alpha = a)
res_MLE


res2_MLE <- results(Deseq2_MLE, name="specific_HighAL_vs_control", alpha = a)
res2_MLE

res3_MLE <- results(Deseq3_MLE, name="specific_HighAL_vs_control", alpha = a)
res3_MLE

summary(res_MLE)
summary(res2_MLE)
summary(res3_MLE)

#Shrinkage

# res2_MLE_shrink <- lfcShrink(Deseq2_MLE, coef=2, type = "apeglm")
# 
# plotMA(res2_MLE)
# plotMA(res2_MLE_shrink)

# LFC shrinkage 
# Additional step to enhance the accuracy of analysis
# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05b_wald_test_results.html

write.csv(as.matrix(res_MLE), file = "res_MLE_DESeq2.csv", row.names = T)
write.csv(as.matrix(res2_MLE), file = "res2_MLE_DESeq2.csv", row.names = T)
write.csv(as.matrix(res3_MLE), file = "res3_MLE_DESeq2.csv", row.names = T)

# manually correct the first column in excel


#?rownames_to_column
as.data.frame(res_MLE) %>% 
  filter(padj<a) -> res_sig_MLE

as.data.frame(res2_MLE) %>% 
  filter(padj<a) -> res2_sig_MLE

as.data.frame(res3_MLE) %>% 
  filter(padj<a) -> res3_sig_MLE


summary(res_sig_MLE)
summary(res2_sig_MLE)
summary(res3_sig_MLE)

# plot of significant features
# res2_sig_MLE_shrink <- lfcShrink(Deseq2_sig_MLE, coef=2, type = "apeglm")
# 
# plotMA(res2_MLE_shrink)

write.csv(as.matrix(res_sig_MLE), file = "sig_res_MLE_DESeq2.csv", row.names = T)
write.csv(as.matrix(res2_sig_MLE), file = "sig_res2_MLE_DESeq2.csv", row.names = T)
write.csv(as.matrix(res3_sig_MLE), file = "sig_res3_MLE_DESeq2.csv", row.names = T)


