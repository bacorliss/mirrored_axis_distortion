

#https://biostatsquid.com/volcano-plots-r-tutorial/

#https://sdgamboa.github.io/post/2020_volcano/

# BiocManager::install("LTLA/scRNAseq")
# install.packages("remotes")




sce <- BacherTCellData(filtered = TRUE, ensembl = FALSE, location = TRUE)


counts <- assay(sce, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)
assayNames(sce)


# Loading relevant libraries 
library(tidyverse) 
library(RColorBrewer) 
library(cowplot) 
library(ggrepel) 
source("R/mirrored_axis_distortion.R")

# Path initialization
base_dir = "mirror_fc"
fig_num = "9" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
data_path <- paste0(base_dir,"/data/")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2.5,2.5)



# Import DGE results
df <- read.csv(paste0(getwd(),"/",data_path, 'severevshealthy_degresults.csv'), row.names = 1)


# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$log2fc > 0.6 & df$pval < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$log2fc < -0.6 & df$pval < 0.05] <- "DOWN"

df$fc <- 2^df$log2fc
df$mfc <- fc_to_mfc(df$fc)
df$con_mfc <- contract1(df$mfc)


g1 <- ggplot(data = df, aes(x = log2fc, y = -log10(pval), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values = c("blue", "black", "red"), labels = c("Down", "NS", "Up")) +
  labs(color = 'Severe', x = expression(log[2]*"FC"), y = expression(-log[10]*"p-value")) +
  theme_classic(base_size = 8) +
  theme(legend.position = "none")
g1
save_plot(paste(fig_path, '/', "A_volcano_logFC.png", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  




# Fold Change Linear Plot
g2 <- ggplot(data=df, aes(x=fc, y=-log10(pval))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  # geom_vline(xintercept = 1, color = "grey") +
  geom_point(aes(col=diffexpressed), size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8)   +
  ylab(expression(-log[10]~(`p-value`))) +
  xlab("FC") +
  theme(legend.position = "none")     
g2
save_plot(paste(fig_path, '/', "B_volcano_FC.png", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# MAD Fold Change Linear Plot
g3 <- ggplot(data=df, aes(x=con_mfc, y=-log10(pval))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  # geom_vline(xintercept = 0, color = "grey") +
  geom_point(aes(col=diffexpressed), size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8) +
  ylab(expression(-log[10]~(`p-value`))) +
  xlab("MAD-FC") +
  theme(legend.position = "none") +
  # scale_x_continuous(breaks=c(c(-9, - 4, 0, seq(4,30,5))), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))
g3
g3 <- gg_revaxis_mfc(g3,'x', num_format = "fraction")
g3
save_plot(paste(fig_path, '/', "C_volcano_MADFC.png", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

