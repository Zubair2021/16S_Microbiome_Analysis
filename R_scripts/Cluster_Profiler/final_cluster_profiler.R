
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


# reading in data from deseq2
df = read.csv("sig_res2_MLE_DESeq2.csv")
head(df)

# we want the log2 fold change 

original_gene_list <- df$log2FoldChange
head(original_gene_list)

# name the vector
names(original_gene_list) <- df$X

# KO_names <- names(original_gene_list)
# head(KO_names)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)

gene_names <- names(gene_list)
head(gene_names)

# KEGG enrichment analysis object
?gseKEGG

# Choose the right organism or try to play with the universe
gsea <- gseKEGG(geneList     = gene_list,
                organism     = "ko",  #organism ko defines all the kegg orthology organisms
                nPerm        = 10000,
                minGSSize    = 3,
                pAdjustMethod = "BH",
                use_internal_data = F,
                by = "fgsea",
                maxGSSize    = 800,
                pvalueCutoff = 0.05,
                keyType       = "kegg",
)

head(gsea, 10)

#?enrichKEGG()

rm(KO_names)
enrich <- enrichKEGG(
  gene_names,
  keyType = "kegg",
  organism = "ko",
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

head(enrich, 10)

# KEGG Module 
enrich_module <- enrichMKEGG(
  gene_names,
  keyType = "kegg",
  organism = "ko",
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
)
head(enrich_module)

# Generate a dotplot
#dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)


dotplot(gsea)
ggsave("./outputs/gsea_dotplot.png", plot = dotplot(gsea), width = 14, height = 8)

dotplot(enrich)
ggsave("./outputs/enrich_dotplot.png", plot = dotplot(enrich), width = 14, height = 8)

dotplot(enrich_module)
ggsave("./outputs/module_dotplot.png", plot = dotplot(enrich_module), width = 14, height = 8)

# The colors for the following are not produced due to the foldchange
cnetplot(enrich,  foldchange=gene_list)
cnetplot(enrich_module,categorySize="p.adjust", foldchange=original_gene_list)
head(gene_list)

# categorySize="pvalue",

?cnetplot
head(gene_list)


?cnetplot

# Ridgeplot
library(ggridges)
ridgeplot(gsea) 
ggsave("./outputs/gsea_ridgeplot.png", plot = ridgeplot(gsea) , width = 14, height = 8)

# did not work despite trying multiple times

# GSEA plot to visualize the GSEA result
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gsea, by = "all", title = gsea$Description[2], geneSetID = 1)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=gene_list, pathway.id="02020", species = "ko")
dme

pathview(gene.data=gene_list, pathway.id="02020", species = "ko")

# Produce a different plot (PDF) (not displayed here)
pathview(gene.data=gene_list, pathway.id="01230", species = "ko", kegg.native = F)

#install.packages("ggupset")
library("ggupset")

png(file = "./output/upset_plot.png")
upsetplot(gsea)
dev.off()

ggsave("./outputs/gsea_upsetplot.png", plot = upsetplot(gsea) , width = 14, height = 8)


# barplot 
obj <- barplot(enrich, 
        drop = T, 
        showCategory = 20, 
        title = "Significantly Regulated Pathways",
        font.size = 12)

ggsave("./outputs/enrich_barplot.png", obj, width = 14, height = 8)

dev.off()


# need to work on it before I proceed
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/
# Enrich plot

