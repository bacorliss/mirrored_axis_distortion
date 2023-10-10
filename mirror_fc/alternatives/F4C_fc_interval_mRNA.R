
# https://arthritis-research.biomedcentral.com/articles/10.1186/ar2053
#https://arthritis-research.biomedcentral.com/articles/10.1186/ar2053/figures/1

# Mislabeled axis



library(tidyverse)
library(cowplot)
source("R/mirrored_axis_distortion.R")

base_dir = "mirror_fc"
fig_num = "4" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2.5,2.25)



df0 <- read.csv(paste0(base_dir, "/mrna_expression_cartilege_Clements2006_F1.csv"))
df0$x <- factor(df0$x, order = TRUE)
df <- df0
  
# Although research paper labels axis as fold change, its actually relative change
# Need to convert relative change data to fold change
df[,c("hi","mid","lo")] <- df0[,c("hi","mid","lo")] + 1

df$log2_hi <- log2(df$hi)
df$log2_mid <- log2(df$mid)
df$log2_lo <- log2(df$lo)

df$mfc_hi <- contract1(fc_to_mfc(df$hi))
df$mfc_mid <- contract1(fc_to_mfc(df$mid))
df$mfc_lo <- contract1(fc_to_mfc(df$lo))


# Log2 Plot
g1 <- ggplot(data = df,aes(x=x, y = log2_mid)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(ymax = log2_hi, ymin = log2_lo), linewidth = .6, width = .6) + 
  geom_point(size = 1.1, alpha = 1) +
  theme_classic(base_size = 8) + 
  scale_x_discrete(labels = df$label) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7),
        legend.position = "none"
  ) +
  scale_y_continuous(oob=squish) +
  # coord_cartesian(ylim = c(-4,2)) +
  xlab("") + ylab("Log2(FC) From Vehicle Tx")
g1
save_plot(paste(fig_path, '/', "A_gene-int_log2.png", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# FC Plot
g2 <- ggplot(data = df,aes(x=x, y = mid)) +
  geom_hline(yintercept = 1, color = "grey") +
  geom_errorbar(aes(ymax = hi, ymin = lo), linewidth = .6, width = .6) + 
  geom_point(size = 1.1,alpha = 1) +
  theme_classic(base_size = 8) + 
  scale_x_discrete(labels = df$label) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7),
        legend.position = "none") +
  # scale_y_continuous(expand = c(0,0),breaks=seq(0,3,.5)) +
  xlab("") + ylab("FC From Vehicle Tx")
g2
save_plot(paste(fig_path, '/', "B_gene-int_fc.png", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# MAD FC Plot
g3 <- ggplot(data = df,aes(x=x, y = mfc_mid)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(ymax = mfc_hi, ymin = mfc_lo), linewidth = .6, width = .6) + 
  geom_point(size = 1.2,alpha = 1) +
  theme_classic(base_size = 8) + 
  scale_x_discrete(labels = df$label) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7),
        legend.position = "none") +
  # scale_y_continuous(expand = c(0,0)) +
  # coord_cartesian(ylim = c(-9,3)) +
  xlab("") + ylab("MAD-FC From Vehicle Tx")
g3
g3 <- gg_revaxis_mfc(g3,'y', num_format = "fraction")
g3
save_plot(paste(fig_path, '/', "c_gene-int_mad-fc.png", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  