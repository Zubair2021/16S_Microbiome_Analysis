# Install & load dependencies

#install.packages(c('lmerTest','pbapply','car','dplyr','vegan','chemometrics','ggplot2','pheatmap','hash','logging','data.table','MuMIn','glmmTMB','MASS','cplm','pscl'), repos='http://cran.r-project.org')

dependencies <- c('lmerTest','pbapply','car','dplyr','vegan','chemometrics','ggplot2','pheatmap','hash','logging','data.table','MuMIn','glmmTMB','MASS','cplm','pscl')
lapply(dependencies, library, character.only=TRUE)

# Package installation

# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Maaslin2")

library(Maaslin2)
input_data <- read.csv("transposed-genus-table.csv")
head(input_data)



# Example run

?Maaslin2

maaslin2_fit <- Maaslin2(
  input_data = "transposed-genus-table.txt",
  input_metadata = "metadata.txt",
  "histo_new",
  min_abundance = 0,
  min_prevalence = 0.1,
  min_variance = 0.0,
  max_significance = 0.05,
  normalization = 'TSS',
  transform = "AST", # because it is also sensitive to zeros, log is not
  analysis_method = "LM",
  fixed_effects = c('Treatment', 'Days_Post_inoculation', 'Body_Weight', 'Crypt_Depth', 'Crypt_Width'),
  random_effects = 'ID',
  correction = "BH",
  standardize = TRUE,
  cores = 8,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = 'Treatment, Control'
)
