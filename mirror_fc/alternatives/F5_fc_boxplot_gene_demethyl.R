

# https://www.frontiersin.org/articles/10.3389/fimmu.2017.00412/full
# Fig 4
# Relationship between DMRs and fold change of genes. Boxplot depicting log2-fold
# change in differentially expressed genes (≥4-fold change) and their DMRs in each 
# and whole gene region. A group of randomly selected genes is also shown (random)



library(tidyverse)
library(cowplot)
source("R/mirrored_axis_distortion.R")

base_dir = "mirror_fc"
fig_num = "5" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2.5,2.25)



df0 <- read.csv(paste0(base_dir, "/boxplot_genebody_demethylation_zou_2020.csv"))
df_log2 <- df0
df_log2$x <- factor(df0$x, order = TRUE)
df_log2$group <- factor(df0$group, order = TRUE)
df_log2$dir <- factor(df0$dir, order = TRUE)


df_lin <- df_log2
df_lin[,2:6] <- sapply(df_log2[,2:6], function(x) 2^x)


df_mad <- df_lin
df_mad[,2:6]<-  sapply(df_lin[,2:6], function(x) contract1(fc_to_mfc(x))) 


# Although research paper labels axis as fold change, its actually relative change
# Need to convert relative change data to fold change





##  Log Plot
g0 <- ggplot(df_log2, aes(x=x, ymin = min, lower = lower, middle = mid, 
                          upper = upper, ymax = max,
                          fill = dir, color = x)) +
  geom_hline(yintercept = 0, size = .3, color = "grey70") +
  geom_boxplot(stat = "identity", position = "identity") +
  # scale_y_continuous(breaks = seq(-3,3,1),
  #                    labels = seq(-3,3,1),limits = c(-3.5,3.5)) + 
  # scale_x_continuous(breaks = df_ints$fcu,
  #                    labels = df_ints$fcu) +
  scale_fill_manual(values = rep("white",5)) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab(expression(log[2]~FC)) + xlab("FC from No Change") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g0
save_plot(paste(fig_path, '/', "A_boxplot_log2-fc.jpg", sep = ""),
          g0, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  







##  Linear Plot
g1 <- ggplot(df_lin, aes(x=x, ymin = min, lower = lower, middle = mid, 
                         upper = upper, ymax = max,
                         fill = dir, color = x)) +
  geom_hline(yintercept = 0, size = .3, color = "grey70") +
  geom_boxplot(stat = "identity", position = "identity") +
  # scale_y_continuous(breaks = seq(-3,3,1),
  #                    labels = seq(-3,3,1),limits = c(-3.5,3.5)) + 
  # scale_x_continuous(breaks = df_ints$fcu,
  #                    labels = df_ints$fcu) +
  scale_fill_manual(values = rep("white",5)) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab(expression(log[2]~FC)) + xlab("FC from No Change") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g1
save_plot(paste(fig_path, '/', "B_boxplot_lin-fc.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  




##  Log Plot
g2 <- ggplot(df_mad, aes(x=x, ymin = min, lower = lower, middle = mid, 
                         upper = upper, ymax = max,
                         fill = dir, color = x)) +
  geom_hline(yintercept = 0, size = .3, color = "grey70") +
  geom_boxplot(stat = "identity", position = "identity") +
  scale_y_continuous(breaks = seq(-10,10,2),
                     labels = seq(-10,10,2),limits = c(-11.5,7)) +
  # scale_x_continuous(breaks = df_ints$fcu,
  #                    labels = df_ints$fcu) +
  scale_fill_manual(values = rep("white",5)) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab("MAD-FC") + xlab("FC from No Change") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g2
g2 <- gg_revaxis_mfc(g2,'y', num_format = "fraction")
g2
save_plot(paste(fig_path, '/', "C_boxplot_mad-fc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  




