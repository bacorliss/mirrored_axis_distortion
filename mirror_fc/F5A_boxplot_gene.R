



# Code used from this tutorial:
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/77-facilitating-exploratory-data-visualization-application-to-tcga-genomic-data/



# Install required Base Packages
base_packages <- c("ggplot2", "tidyverse", "cowplot","BiocManager","ggpubr")
install.packages(setdiff(base_packages, rownames(installed.packages())))  
# Install required Bioconductor Packages
biocm_packages <- c("RTCGA", "RTCGA.mRNA")
bioc_installs <- setdiff(biocm_packages, rownames(installed.packages()))
if (length(bioc_installs)) {BiocManager::install(bioc_installs) }

# Load base packages
lapply(base_packages, library, character.only = TRUE)
# Load Bioconductor packages packages
lapply(biocm_packages, library, character.only = TRUE)
source("R/mirrored_axis_distortion.R")


# Linear visualization
base_dir = "mirror_fc"
fig_num = "5" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2,2)



# Extract expression values from the following datasets
# 1. Breast invasive carcinoma (BRCA.mRNA) 
# 2. Ovarian serous cystadenocarcinoma (OV.mRNA) 
# 3. Lung squamous cell carcinoma (LUSC.mRNA)

# Expression for 5 genes extracted: GATA3, PTEN, XBP1, ESR1 and MUC1 
expr <- expressionsTCGA(BRCA.mRNA, OV.mRNA, LUSC.mRNA,
                        extract.cols = c("GATA3", "PTEN", "XBP1","ESR1", "MUC1"))
expr

expr$dataset <- gsub(pattern = ".mRNA", replacement = "",  expr$dataset)
expr$bcr_patient_barcode <- paste0(expr$dataset, c(1:590, 1:561, 1:154))
expr



ggsize <- c(2.25,2.25)

# Convert dataframe of RNA expresison from wide to long format
df_exp <- expr %>% pivot_longer(cols = !c(bcr_patient_barcode, dataset), names_to = "Gene", values_to = "log2_fc")
df_exp$fc <- 2^df_exp$log2_fc
df_exp$mad_fc <- contract1(fc_to_mfc(df_exp$fc))



# Log FC Plot
g0 <- ggplot(data = df_exp, aes(x = Gene, y = log2_fc, color = dataset)) +
  geom_hline(yintercept = 0, size = .3, color = "grey70") +
  geom_boxplot(size =.3, outlier.size = .3,outlier.alpha = 0.5) +
  ylab(expression(log[2]~(FC))) + xlab("Gene Name") +
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
  ylab("FC") + xlab("Gene Name") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g2
g2 <- gg_revaxis_mfc(g2,'y', num_format = "fraction")
g2
save_plot(paste(fig_path, '/', "boxplot_TCGA_C_mad-fc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2]) 