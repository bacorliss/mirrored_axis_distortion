

# Code and data appropriated from
# https://bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html

# 5.0 Workflow functions for the entire analysis
# The package contains workflow functions that entail the complete analysis and generate a report.
# 
# 5.1LFQ-based DEP analysis
# Differential enrichment analysis of label-free proteomics data can be performed using the LFQ workflow function.



# Install required Base Packages
base_packages <- c("tidyverse", "cowplot","BiocManager", "reshape2")
install.packages(setdiff(base_packages, rownames(installed.packages())))  
# Install required Bioconductor Packages
biocm_packages <- c("DESeq2", "DEP")
bioc_installs <- setdiff(biocm_packages, rownames(installed.packages()))
if (length(bioc_installs)) {BiocManager::install(bioc_installs) }

# Load base packages
lapply(base_packages, library, character.only = TRUE)
# Load Bioconductor packages packages
lapply(biocm_packages, library, character.only = TRUE)
source("R/mirrored_axis_distortion.R")


# Code Initialization
base_dir = "mirror_fc"
fig_num = "7" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2.25,5)


# The data is provided with the package 
data <- UbiLength
experimental_design <- UbiLength_ExpDesign

# The wrapper function performs the full analysis on dataset
data_results <- LFQ(data, experimental_design, fun = "MinProb", 
                    type = "control", control = "Ctrl", alpha = 0.05, lfc = 1)

# Extract the results table from the LFQ analysis
results_table <- data_results$results


# Extract the sign object
full_data <- data_results$dep
# Used a gene expression heatmap function with specialized algorithm for visualizing most interesting genes.
heatmap_fc <- plot_heatmap(full_data, type = "contrast", kmeans = TRUE,
             k = 6, col_limit = 4, show_row_names = FALSE)
# Extract fold changes from heatmap
def_heatmap <- heatmap_fc@ht_list$`log2 Fold change`@matrix


# Extracted log fold changes from results
df_expr <- results_table %>% select(ends_with("_ratio"))
# Added gene names as seprate column
df_expr$Gene <- results_table$name
# Removed extra characteris is sample (column) names
colnames(df_expr) <- str_replace(colnames(df_expr),"_vs_Ctrl_ratio","")
# Discarded genes not included in default heatmap
df_sub <- df_expr[df_expr$Gene %in% rownames(def_heatmap),]
# COnverted genes to factor to ensure same gene order as in default heatmap
df_sub$Gene <- factor(df_sub$Gene, ordered = TRUE, levels = rownames(def_heatmap))
# COnverted wide dataset to long
df_sub <- pivot_longer(df_sub,cols = -Gene, names_to = "Sample", values_to = "log2FC")

# Added columsn for FC and madFC
df_sub$FC <- 2^ df_sub$log2FC
df_sub$madFC <- contract1(fc_to_mfc(df_sub$FC))


ggsize <- c(5,2.5)
saturate <- function(df, colname,lo, hi) {
  df[[colname]][df[[colname]]<lo] <- lo
  df[[colname]][df[[colname]]>hi] <- hi
  return(df)
}

# Only label 1 in 10 genes
gene_list <- rownames(def_heatmap)
gene_sublist <- rep(" ",length(gene_list))
gene_sublist[seq(1,length(gene_list), 10)] <- gene_list[seq(1,length(gene_list), 10)]


# Log2 FC plot
g1 <- ggplot(data = saturate(df_sub,'log2FC',-4,4), aes(x=Sample, y= Gene, fill = log2FC)) + 
  geom_tile() +
  scale_y_discrete("", labels = gene_sublist) +
  scale_fill_gradientn(name = expression(log[2](FC)),limits = c(-4,4),
                       colors=c("#377eb8", "white", "#e41a1c"),
                       # breaks = c(-15,-7, 0, 7, 15), labels = mad_fc_labeller,
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                              barwidth = grid::unit(.6, "npc"),
                                              barheight = grid::unit(.025, "npc")),) +
  theme_classic(base_size = 8) +
  theme(legend.position="top")
g1
save_plot(paste(fig_path, '/', "heatmap_protein_A_log2-fc.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2]) 

# Linear FC plot
g2 <- ggplot(data = saturate(df_sub,'FC',1/16,16), aes(x=Sample, y= Gene, fill = FC)) + 
  geom_tile() +
  scale_y_discrete("", labels = gene_sublist) + 
  scale_fill_gradientn(name = "FC", limits = c(0,16),
                       colors=c("#377eb8", "white", "#e41a1c"), values=c(0, .75/16, 1.25/16, 1),
                       breaks = c(0,-1/16,1, 4,8,12,16),
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                              barwidth = grid::unit(.6, "npc"),
                                              barheight = grid::unit(.025, "npc")),) +
  theme_classic(base_size = 8) +
  theme(legend.position="top")
g2
save_plot(paste(fig_path, '/', "heatmap_protein_B_fc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

# MAD-FC plot
g3 <- ggplot(data = saturate(df_sub,'madFC',-16,16), aes(x=Sample, y= Gene, fill = madFC)) + 
  geom_tile() +
  scale_y_discrete("", labels = gene_sublist) + 
  scale_fill_gradientn(name = "FC", limits = c(-16,16), 
                       colors=c("#377eb8", "white", "#e41a1c"),
                       breaks = c(-15,-7, 0, 7, 15), labels = mad_fc_labeller,
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                               barwidth = grid::unit(.6, "npc"), 
                                              barheight = grid::unit(.025, "npc")),) +
  theme_classic(base_size = 8) +
  theme(legend.position="top")
g3
save_plot(paste(fig_path, '/', "heatmap_protein_C_mad-fc.jpg", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


