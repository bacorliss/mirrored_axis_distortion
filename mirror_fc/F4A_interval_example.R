

library(ggplot2)
library(cowplot)
source("R/mirrored_axis_distortion.R")

# Linear visualization
base_dir = "mirror_fc"
fig_num = "4" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2,2)

df_ints <- data.frame(fcu = seq(-8,8,2), 
                      mid = c(1/9, 1/7, 1/5, 1/3, seq(1,9,2)))
df_ints$lo <- c(1/11, df_ints$mid[1:(nrow(df_ints)-1)])
df_ints$hi <- c(df_ints$mid[2:nrow(df_ints)], 11)

df_ints$log2_mid <- log2(df_ints$mid) 
df_ints$log2_lo  <- log2(df_ints$lo) 
df_ints$log2_hi  <- log2(df_ints$hi) 

df_ints$mad_mid <- contract1(fc_to_mfc(df_ints$mid))
df_ints$mad_lo  <- contract1(fc_to_mfc(df_ints$lo))
df_ints$mad_hi  <- contract1(fc_to_mfc(df_ints$hi))



##  Log Plot
g0 <- ggplot( data = df_ints, aes(x = fcu, y = log2_mid, ymax = log2_hi, ymin = log2_lo,color = fcu)) +
  geom_hline(yintercept = 0, linewidth = .3, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = .3, color = "grey70") +
  geom_errorbar(width = 0.5, size = .4) +
  geom_point(size = .4) + 
  scale_colour_gradient(low = "#e41a1c", high = "#377eb8", guide = "colourbar",
                        aesthetics = "colour") +
  scale_y_continuous(breaks = seq(-3,3,1),
                     labels = seq(-3,3,1),limits = c(-3.5,3.5)) + 
  scale_x_continuous(breaks = df_ints$fcu,
                     labels = df_ints$fcu) +
  ylab(expression(log[2]~(FC))) + xlab("Fold Change Units") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g0
save_plot(paste(fig_path, '/', "A_int_est_log2-fc.jpg", sep = ""),
          g0, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


##  Linear Plot
g1 <- ggplot( data = df_ints, aes(x = fcu, y = mid, ymax = hi, ymin = lo, color = fcu)) +
  geom_hline(yintercept = 1, color = "grey70", linewidth = .3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = .3) +
  geom_errorbar(width = 0.5, size = .3) +
  geom_point(size = .4) + 
  scale_colour_gradient(low = "#e41a1c", high = "#377eb8", guide = "colourbar",
                        aesthetics = "colour") +
  scale_y_continuous(limits = c(0,11)) +
  scale_x_continuous(breaks = df_ints$fcu,
                     labels = df_ints$fcu) +
  ylab("FC") + xlab("Fold Change Units") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g1
save_plot(paste(fig_path, '/', "B_int_est_linear-fc.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


##  MAD Plot
g2 <- ggplot( data = df_ints, aes(x = fcu, y = mad_mid, ymax = mad_hi, ymin = mad_lo, color = fcu)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = .3) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = .3) +
  geom_errorbar(width = 0.5, size = .3) +
  geom_point(size = .4) +
  scale_colour_gradient(low = "#e41a1c", high = "#377eb8", guide = "colourbar",
                        aesthetics = "colour") +
  scale_y_continuous(breaks = seq(-10,10,2),
                     labels = seq(-10,10,2), limits = c(-10.5,10.5)) +
  scale_x_continuous(breaks = df_ints$fcu,
                     labels = df_ints$fcu) +
  ylab("FC") + xlab("Fold Change Units") +
  theme_classic(base_size = 8) + 
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g2
g2 <- gg_revaxis_mfc(g2,'y', num_format = "fraction")
g2
save_plot(paste(fig_path, '/', "C_int_est_mad-fc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

