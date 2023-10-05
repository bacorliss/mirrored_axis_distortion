# F3_FoldChange_Protein_expression

#https://www.oncotarget.com/article/19461/text/

library(tidyverse)
library(cowplot)
source("R/mirrored_axis_distortion.R")

base_dir = "mirror_fc"
fig_num = "4" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(3,2.25)



df0 <- read.csv(paste0(base_dir, "/protein_expression_John O’Shea.csv"))
df <- pivot_longer(df0,cols = -c("X", "label"), 
                   names_to = "color_pt",values_to = "fc") %>% 
  separate(color_pt, c("color", "pt"), '_') %>%
  pivot_wider(names_from = pt, values_from = fc)

df$log2_hi <- log2(df$hi)
df$log2_mid <- log2(df$mid)
df$log2_lo <- log2(df$lo)

df$mfc_hi <- contract1(fc_to_mfc(df$hi))
df$mfc_mid <- contract1(fc_to_mfc(df$mid))
df$mfc_lo <- contract1(fc_to_mfc(df$lo))


# Log2 Plot
g1 <- ggplot(data = df,aes(x=label, y = log2_mid, color = color)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_hline(yintercept = c(-.26, 0.26), color = "grey", linetype = "dashed") +
  geom_errorbar(aes(ymax = log2_hi, ymin = log2_lo), linewidth = .8, width = .8, 
                position=position_dodge(width = 0.6), alpha = 1) + 
  geom_point(size = .8,position=position_dodge(width = 0.6),alpha = 1) +
  theme_classic(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7),
        legend.position = "none"
  ) +
  scale_color_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c("#1b9e77", "#d95f02", "#7570b3")) +
  scale_shape_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c(18, 15, 17)) +
  # scale_y_continuous(expand = c(0,0)) +
  # coord_cartesian(ylim = c(-4,2)) +
  xlab("") + ylab("Log2(FC) From Vehicle Tx")
g1
save_plot(paste(fig_path, '/', "A_gene-int_log2.png", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# FC Plot
g2 <- ggplot(data = df,aes(x=label, y = mid, color = color)) +
  geom_hline(yintercept = 1, color = "grey") +
  geom_hline(yintercept = c(1.2, 0.8), color = "grey", linetype = "dashed") +
  geom_errorbar(aes(ymax = hi, ymin = lo), linewidth = .8, width = .8, 
                position=position_dodge(width = 0.6), alpha = 1) + 
  geom_point(size = .8,position=position_dodge(width = 0.6),alpha = 1) +
  theme_classic(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7),
          legend.position = "none",
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6)
        ) +
  scale_color_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), 
                     values = c("#1b9e77", "#d95f02", "#7570b3")) +
  scale_y_continuous(breaks=seq(0,3,.5)) +
  xlab("") + ylab("FC From Vehicle Tx")
g2
save_plot(paste(fig_path, '/', "B_gene-int_fc.png", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# MAD FC Plot
g3 <- ggplot(data = df,aes(x=label, y = mfc_mid, color = color)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_hline(yintercept = c(.2, -0.2), color = "grey", linetype = "dashed") +
  geom_errorbar(aes(ymax = mfc_hi, ymin = mfc_lo), linewidth = .8, width = .8, 
                position=position_dodge(width = 0.6), alpha = 1) + 
  geom_point(size = .8,position=position_dodge(width = 0.6),alpha = 1) +
  theme_classic(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7),
        legend.position = "none",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
  ) +
  scale_color_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c("#1b9e77", "#d95f02", "#7570b3")) +
  scale_shape_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c(18, 15, 17)) +
  # scale_y_continuous(expand = c(0,0)) +
  # coord_cartesian(ylim = c(-9,3)) +
  xlab("") + ylab("MAD-FC From Vehicle Tx")
g3
g3 <- gg_revaxis_mfc(g3,'y', num_format = "fraction")
g3

save_plot(paste(fig_path, '/', "c_gene-int_mad-fc.png", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  