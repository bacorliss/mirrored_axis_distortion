


# Code appropriated from
# https://sdgamboa.github.io/post/2020_volcano/
# Processed data found at:
# https://github.com/sdgamboa/misc_datasets
# Data originally collected from
#  S. D. Gamboa-Tuz, A. Pereira-Santana, J. A. Zamora-Briseño, E. Castano, F. Espadas-Gil, J. T. Ayala-Sumuano, M. Á. Keb-Llanes, F. Sanchez-Teyer, L. C. Rodríguez-Zapata, Transcriptomics and co-expression networks reveal tissue-specific responses and regulatory hubs under mild and severe drought in papaya (Carica papaya L.). Sci Rep. 8, 14539 (2018).

# BiocManager::install("LTLA/scRNAseq")
# install.packages("remotes")



library(tidyverse)
library(ggrepel)
library(cowplot)
source("R/mirrored_axis_distortion.R")


# Path initialization
base_dir = "mirror_fc"
fig_num = "8" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
data_path <- paste0(base_dir,"/data/")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2.5,2.5)



# Import data
df <- read_tsv("https://raw.githubusercontent.com/sdgamboa/misc_datasets/master/L0_vs_L20.tsv")


# Annotated up and down regulated
df <- df %>%  
  mutate(Expression = case_when(logFC >= 1 & FDR <= 0.05 ~ "Up",
                           logFC <= -1 & FDR <= 0.05 ~ "Down",
                           TRUE ~ "NS"))


# Calculate FOld change, Contracted fold change, MAD-FC
df$fc <- 2^df$logFC
df$mfc <- fc_to_mfc(df$fc)
df$con_mfc <- contract1(df$mfc)








g1 <- ggplot(data = df, aes(x = logFC, y = -log(FDR,10), col = Expression)) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values = c("blue", "black", "red"), labels = c("Down", "NS", "Up")) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab(expression(-log[10]~(`FDR`))) +
  xlab(expression(log[2]*("FC"))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none")
g1
save_plot(paste(fig_path, '/', "A_volcano_logFC.png", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

# Fold Change Linear Plot
g2 <- ggplot(data = df, aes(x = fc, y = -log(FDR,10), col = Expression)) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values = c("blue", "black", "red"), labels = c("Down", "NS", "Up")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 8)   +
  ylab(expression(-log[10]~(`FDR`))) +
  xlab("FC") +
  theme(legend.position = "none")     
g2
save_plot(paste(fig_path, '/', "B_volcano_FC.png", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# MAD Fold Change Linear Plot
g3 <- ggplot(data=df, aes(x=con_mfc, y= -log(FDR,10), col = Expression)) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values=c("blue", "black", "red"), labels = c("Down", "NS", "Up")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 8) +
  ylab(expression(-log[10]~(`FDR`))) +
  xlab("FC") +
  theme(legend.position = "none") +
  # scale_x_continuous(breaks=c(c(-9, - 4, 0, seq(4,30,5))), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))
g3
g3 <- gg_revaxis_mfc(g3,'x', num_format = "fraction")
g3
save_plot(paste(fig_path, '/', "C_volcano_MADFC.png", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  
















