

library(ggplot2)
library(cowplot)
source("R/mirrored_axis_distortion.R")

# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/77-facilitating-exploratory-data-visualization-application-to-tcga-genomic-data/




# Linear visualization
base_dir = "mirror_fc"
fig_num = "5" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2,2)


# install.packages("ggpubr")
library(ggpubr)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")


# Install the main RTCGA package

library(BiocManager)
# BiocManager::install("RTCGA")
# # Install the clinical and mRNA gene expression data packages
# BiocManager::install("RTCGA.clinical")
# BiocManager::install("RTCGA.mRNA")



# Preprocessing
library(RTCGA)
library(RTCGA.mRNA)
expr <- expressionsTCGA(BRCA.mRNA, OV.mRNA, LUSC.mRNA,
                        extract.cols = c("GATA3", "PTEN", "XBP1","ESR1", "MUC1"))
expr

expr$dataset <- gsub(pattern = ".mRNA", replacement = "",  expr$dataset)
expr$bcr_patient_barcode <- paste0(expr$dataset, c(1:590, 1:561, 1:154))
expr


library(ggplot2)
library(tidyverse)
library(cowplot)
source("R/mirrored_axis_distortion.R")

ggsize <- c(2.25,2.25)

df_exp <- expr %>% pivot_longer(cols = !c(bcr_patient_barcode, dataset), names_to = "Gene", values_to = "log2_fc")
df_exp$fc <- 2^df_exp$log2_fc
df_exp$mad_fc <- contract1(fc_to_mfc(df_exp$fc))



# Log FC Plot
g0 <- ggplot(data = df_exp, aes(x = Gene, y = log2_fc, color = dataset)) +
  geom_hline(yintercept = 0, size = .3, color = "grey70") +
  geom_boxplot(size =.3, outlier.size = .3,outlier.alpha = 0.5) +
  ylab(expression(log[2]~FC)) + xlab("Gene Name") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g0
save_plot(paste(fig_path, '/', "boxplot_TCGA_A_log2-fc.jpg", sep = ""),
          g0, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

# Linear FC Plot
g1 <- ggplot(data = df_exp, aes(x = Gene, y = fc, color = dataset)) +
  geom_hline(yintercept = 1, size = .3, color = "grey70") +
  geom_boxplot(size =.3, outlier.size = .3,outlier.alpha = 0.5) +
  ylab("FC") + xlab("Gene Name") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g1
save_plot(paste(fig_path, '/', "boxplot_TCGA_B_fc.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

# MAD FC Plot
g2 <- ggplot(data = df_exp, aes(x = Gene, y = mad_fc, color = dataset)) +
  geom_hline(yintercept = 1, size = .3, color = "grey70") +
  geom_boxplot(size =.3, outlier.size = .3,outlier.alpha = 0.5) +
  scale_y_continuous(breaks = seq(-70,70,10),
                     labels = seq(-70,70,10)) +
  ylab("MAD-FC") + xlab("Gene Name") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g2
g2 <- gg_revaxis_mfc(g2,'y', num_format = "fraction")
g2
save_plot(paste(fig_path, '/', "boxplot_TCGA_C_mad-fc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  




# Violin plots of same data
fig_num = "6" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)



# Log Plot
g0 <- ggplot(data = df_exp, aes(x = Gene, y = log2_fc, color = dataset)) +
  geom_hline(yintercept = 0, size = .3, color = "grey70") +
  geom_violin(scale = "width", size =.2) +
  ylab(expression(log[2]~FC)) + xlab("Gene Name") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g0
save_plot(paste(fig_path, '/', "violin_TCGA_A_log2-fc.jpg", sep = ""),
          g0, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

# Linear Plot
g1 <- ggplot(data = df_exp, aes(x = Gene, y = fc, color = dataset)) +
  geom_hline(yintercept = 1, color = "grey70") +
  geom_violin(scale = "width", size = .2) +
  ylab("FC") + xlab("Gene Name") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g1
save_plot(paste(fig_path, '/', "violin_TCGA_B_fc.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

# MAD FC Plot
g2 <- ggplot(data = df_exp, aes(x = Gene, y = mad_fc, color = dataset)) +
  geom_hline(yintercept = 1, color = "grey70") +
  geom_violin(scale = "width", size = .2) +
  scale_y_continuous(breaks = seq(-60,60,20),
                     labels = seq(-60,60,20)) +
  ylab("MAD-FC") + xlab("Gene Name") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g2
g2 <- gg_revaxis_mfc(g2,'y', num_format = "fraction")
g2
save_plot(paste(fig_path, '/', "violin_TCGA_C_mad-fc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  



