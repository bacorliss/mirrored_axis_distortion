

library(ggplot2)
library(cowplot)
source("R/mirrored_axis_distortion.R")

# Linear visualization
base_dir = "mirror_fc"
fig_num = "6" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2,2)



mad_fc_gen <- rnorm(1000, mean = 0, sd = .5)
mad_fc_gen <- mad_fc_gen * (1+ abs(mad_fc_gen))^.5

fcu <- seq(-8,8,2)

df_list = list()
for (x in seq_along(fcu)) {
  df_list[[x]] <- data.frame(x = fcu[x], mad_fc = mad_fc_gen+fcu[x])
}

df_fc <- do.call(rbind, df_list)
df_fc$fc <- mfc_to_fc(rev_contract1(df_fc$mad_fc))
df_fc$log2_fc <- log2(df_fc$fc)
# df_fc$x <- factor(df_fc$x, ordered =  TRUE)


##  Log FC Plot, violin
g0 <- ggplot(df_fc, aes(x=factor(x), y=log2_fc, color = x)) +
  geom_hline(yintercept = 0, size = .3, color = "grey70") +
  geom_violin(scale = "width", size =.2) +
  scale_colour_gradient(low = "#e41a1c", high = "#377eb8", guide = "colourbar",
                        aesthetics = "colour") +
  ylab(expression(log[2]~(FC))) + xlab("Fold Change Units") +
  theme_classic(base_size = 8) +
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g0
save_plot(paste(fig_path, '/', "A_violin_example_log2-fc.jpg", sep = ""),
          g0, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


##  Log FC Plot, violin
g1 <- ggplot(df_fc, aes(x=factor(x), y=fc, color = x)) +
  geom_hline(yintercept = 1, size = .3, color = "grey70") +
  geom_violin(scale = "width", size =.2) +
  scale_colour_gradient(low = "#e41a1c", high = "#377eb8", guide = "colourbar",
                        aesthetics = "colour") +
  ylab("FC") + xlab("Fold Change Units") +
  theme_classic(base_size = 8) +
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g1
save_plot(paste(fig_path, '/', "B_violin_example_fc.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


## MAD-FC Plot, violin
g2 <- ggplot(df_fc, aes(x=factor(x), y=mad_fc, color = x)) +
  geom_hline(yintercept = 1, size = .3, color = "grey70") +
  geom_violin(scale = "width", size =.2) +
  scale_colour_gradient(low = "#e41a1c", high = "#377eb8", guide = "colourbar",
                        aesthetics = "colour") +
  ylab("FC") + xlab("Fold Change Units") +
  theme_classic(base_size = 8) +
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g2
save_plot(paste(fig_path, '/', "C_violin_example_mad-fc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

