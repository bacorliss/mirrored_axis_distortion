

# https://bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html

# 5.0 Workflow functions for the entire analysis
# The package contains workflow functions that entail the complete analysis and generate a report.
# 
# 5.1LFQ-based DEP analysis
# Differential enrichment analysis of label-free proteomics data can be performed using the LFQ workflow function.


library(DESeq2)
library(ggplot2)
library(DEP)
source("R/mirrored_axis_distortion.R")

# The data is provided with the package 
data <- UbiLength
experimental_design <- UbiLength_ExpDesign

# The wrapper function performs the full analysis
data_results <- LFQ(data, experimental_design, fun = "MinProb", 
                    type = "control", control = "Ctrl", alpha = 0.05, lfc = 1)



# Extract the results table
results_table <- data_results$results

# Number of significant proteins
results_table %>% filter(significant) %>% nrow()


# Extract the sign object
full_data <- data_results$dep

# Use the full data to generate a heatmap
heatmap_fc <-plot_heatmap(full_data, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE)

mx_log2_fc <- heatmap_fc@ht_list$`log2 Fold change`@matrix
mx_fc <- 2^mx_log2_fc

mx_mad_fc <- mx_fc
mx_mad_fc[] <- vapply(mx_mad_fc, function(x) contract1(fc_to_mfc(x)), numeric(1))

library(tidyverse)
library(reshape2)


df_log2 <- as.data.frame(melt(mx_log2_fc, varnames=c('Gene', 'Sample'), value.name = "y"))
df_log2$Gene <- factor(df_log2$Gene, ordered = TRUE, levels = rownames(mx_log2_fc))
levels(df_log2$Sample) <- c("Ubi1", "Ubi4", "ubi6")

df_lin <- as.data.frame(melt(mx_fc, varnames=c('Gene', 'Sample'), value.name = "y"))
df_lin$Gene <- factor(df_lin$Gene, ordered = TRUE, levels = rownames(mx_log2_fc))
levels(df_lin$Sample) <- c("Ubi1", "Ubi4", "ubi6")

df_mad <- as.data.frame(melt(mx_mad_fc, varnames=c('Gene', 'Sample'), value.name = "y"))
df_mad$Gene <- factor(df_mad$Gene, ordered = TRUE, levels = rownames(mx_log2_fc))
levels(df_mad$Sample) <- c("Ubi1", "Ubi4", "ubi6")

# heatmap.2(mx_log2_fc, dendrogram = "both",scale = "none", trace="none", 
#           density.info="none", col="heat.colors")



saturate <- function(df, colname,lo, hi) {
  df[[colname]][df[[colname]]<lo] <- lo
  df[[colname]][df[[colname]]>hi] <- hi
  return(df)
}

ggplot(data = saturate(df_log2,'y',-4,4), aes(x=Sample, y= Gene, fill = y)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(-4,4), 
                       colours=c("#377eb8", "white", "#e41a1c"))



ggplot(data = saturate(df_lin,'y',1/16,16), aes(x=Sample, y= Gene, fill = y)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(0,16), 
                       colours=c("#377eb8", "white", "#e41a1c"))



ggplot(data = saturate(df_mad,'y',-16,16), aes(x=Sample, y= Gene, fill = y)) + 
  geom_tile() +
  scale_fill_gradientn(limits = c(-16,16), 
                       colours=c("#377eb8", "white", "#e41a1c"))


