

library(ggplot2)
library(cowplot)
source("R/mirrored_axis_distortion.R")

# Linear visualization
base_dir = "mirror_fc"
fig_num = "5" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2,2)

# x       max     upper       mid     lower       min  group  dir
df_lin <- data.frame(x = seq(-8,8,2), 
                      mid = c(1/9, 1/7, 1/5, 1/3, seq(1,9,2)))
df_lin$lower <- c(1/11, df_lin$mid[1:(nrow(df_lin)-1)])
df_lin$min <- c(1/13, df_lin$lower[1:(nrow(df_lin)-1)])
df_lin$upper <- c(df_lin$mid[2:nrow(df_lin)], 11)
df_lin$max <- c(df_lin$upper[2:nrow(df_lin)], 13)


df_log2 <- df_lin
df_log2[,2:6] <- sapply(df_lin[,2:6], function(x) log2(x))


df_mad <- df_lin
df_mad[,2:6] <-  sapply(df_lin[,2:6], function(x) contract1(fc_to_mfc(x))) 



##  Log Plot
g0 <- ggplot(df_log2, aes(x=x, ymin = min, lower = lower, middle = mid, 
                          upper = upper, ymax = max, group = x, color = x)) +
  geom_hline(yintercept = 0, size = .3, color = "grey70") +
  geom_boxplot(stat = "identity", width = 0.75) +
  # scale_y_continuous(breaks = seq(-3,3,1),
  #                    labels = seq(-3,3,1),limits = c(-3.5,3.5)) + 
  # scale_x_continuous(breaks = df_ints$fcu,
  #                    labels = df_ints$fcu) +
  scale_colour_gradient(low = "#e41a1c", high = "#377eb8", guide = "colourbar",
                        aesthetics = "colour") +
  # scale_fill_manual(values = rep("white",5)) +
  # scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab(expression(log[2]~FC)) + xlab("Fold Change Units") +
  theme_classic(base_size = 8) +
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g0
save_plot(paste(fig_path, '/', "A_boxplot_example_log2-fc.jpg", sep = ""),
          g0, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  




#  Log Plot Overlay
# g0b <- ggplot( data = df_log2, aes(x = rep(0,nrow(df_lin)), y = rep(0,nrow(df_lin)),
#                                     ymax = log2_hi-log2_mid, ymin = log2_lo-log2_mid, color = x) ) +
#   # geom_hline(yintercept = 0, size = .3) +
#   geom_vline(xintercept = 0, size = .3, color = "grey70") +
#   geom_boxplot(stat = "identity", width = 0.75) +
#   # scale_y_continuous(breaks = seq(-3,3,1),
#                      # labels = seq(-3,3,1),limits = c(-3.5,3.5)) +
#   # scale_x_continuous(limits = c(-1,1)) +
#   ylab(" ") + xlab(" ") +
#   theme_classic(base_size = 8) +
#   theme(axis.text.x=element_blank(), #remove x axis labels
#               axis.ticks.x=element_blank(), #remove x axis ticks
#               axis.text.y=element_blank(),  #remove y axis labels
#               axis.ticks.y=element_blank(),  #remove y axis ticks
#               legend.position = "none",
#         axis.line.x=element_line(color="white"),
#         axis.line.y=element_line(color="white"))
# g0b
# save_plot(paste(fig_path, '/', "A_int_est_log2-fc_overlay.jpg", sep = ""),
#           g0b, dpi = 600, base_height = ggsize[1], 
#           base_width = .5)  



##  Linear Plot
g1 <- ggplot(df_lin, aes(x=x, ymin = min, lower = lower, middle = mid, 
                          upper = upper, ymax = max, group = x, color = x)) +
  geom_hline(yintercept = 1, size = .3, color = "grey70") +
  geom_boxplot(stat = "identity", width = 0.75) +
  # scale_y_continuous(breaks = seq(-3,3,1),
  #                    labels = seq(-3,3,1),limits = c(-3.5,3.5)) + 
  # scale_x_continuous(breaks = df_ints$fcu,
  #                    labels = df_ints$fcu) +
  scale_colour_gradient(low = "#e41a1c", high = "#377eb8", guide = "colourbar",
                        aesthetics = "colour") +
  # scale_fill_manual(values = rep("white",5)) +
  # scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab("FC") + xlab("Fold Change Units") +
  theme_classic(base_size = 8) +
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g1
save_plot(paste(fig_path, '/', "B_boxplot_example_lin-fc.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


##  Linear Plot Overlay
# g1b <- ggplot( data = df_lin, aes(x = rep(0,nrow(df_lin)), y = rep(0,nrow(df_lin)),
#                                   ymax = hi-mid, ymin = lo-mid, color = fcu) ) +
#   # geom_hline(yintercept = 0, size = .3) +
#   geom_vline(xintercept = 0, size = .3, color = "grey70") +
#   geom_errorbar(width = 0.7, size = .3) +
#   geom_point(size = .4) + 
#   scale_y_continuous(limits = c(-5,6)) + 
#   scale_x_continuous(limits = c(-1,1)) +
#   ylab(" ") + xlab(" ") +
#   theme_classic(base_size = 8) + 
#   theme(axis.text.x=element_blank(), #remove x axis labels
#         axis.ticks.x=element_blank(), #remove x axis ticks
#         axis.text.y=element_blank(),  #remove y axis labels
#         axis.ticks.y=element_blank(),  #remove y axis ticks
#         legend.position = "none",
#         axis.line.x=element_line(color="white"),
#         axis.line.y=element_line(color="white"))
# g1b
# save_plot(paste(fig_path, '/', "B_int_est_linear-fc_overlay.jpg", sep = ""),
#           g1b, dpi = 600, base_height = ggsize[1], 
#           base_width = .5)  




##  MAD Plot
g2 <- ggplot(df_mad, aes(x=x, ymin = min, lower = lower, middle = mid, 
                         upper = upper, ymax = max, group = x, color = x)) +
  geom_hline(yintercept = 0, size = .3, color = "grey70") +
  geom_boxplot(stat = "identity", width = 0.75) +
  # scale_y_continuous(breaks = seq(-3,3,1),
  #                    labels = seq(-3,3,1),limits = c(-3.5,3.5)) + 
  # scale_x_continuous(breaks = df_ints$fcu,
  #                    labels = df_ints$fcu) +
  scale_colour_gradient(low = "#e41a1c", high = "#377eb8", guide = "colourbar",
                        aesthetics = "colour") +
  # scale_fill_manual(values = rep("white",5)) +
  # scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  ylab("MAD-FC") + xlab("Fold Change Units") +
  theme_classic(base_size = 8) +
  theme(panel.grid.minor = element_blank(),legend.position = "none")
g2
g2 <- gg_revaxis_mfc(g2,'y', num_format = "fraction")
g2
save_plot(paste(fig_path, '/', "C_int_est_mad-fc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


##  Linear Plot Overlay
# g2b <- ggplot( data = df_lin, aes(x = rep(0,nrow(df_lin)), y = rep(0,nrow(df_lin)),
#                                    ymax = mad_hi-mad_mid, ymin = mad_lo-mad_mid, color = fcu) ) +
#   # geom_hline(yintercept = 0, size = .3) +
#   geom_vline(xintercept = 0, size = .3, color = "grey70") +
#   geom_errorbar(width = 0.7, size = .3) +
#   geom_point(size = .4) + 
#   scale_y_continuous(limits = c(-10.5,10.5)) + 
#   scale_x_continuous(limits = c(-1,1)) +
#   ylab(" ") + xlab(" ") +
#   theme_classic(base_size = 8) + 
#   theme(axis.text.x=element_blank(), #remove x axis labels
#         axis.ticks.x=element_blank(), #remove x axis ticks
#         axis.text.y=element_blank(),  #remove y axis labels
#         axis.ticks.y=element_blank(),  #remove y axis ticks
#         legend.position = "none",
#         axis.line.x=element_line(color="white"),
#         axis.line.y=element_line(color="white"))
# g2b
# save_plot(paste(fig_path, '/', "C_int_est_mad-fc_overlay.jpg", sep = ""),
#           g2b, dpi = 600, base_height = ggsize[1], 
#           base_width = .5)  



