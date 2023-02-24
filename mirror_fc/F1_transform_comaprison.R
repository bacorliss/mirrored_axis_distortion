

library(ggplot2)
library(cowplot)
library(scales)
source("R/mirrored_axis_distortion.R")

# Linear visualization
base_dir = "mirror_fc"
fig_num = "1" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

ggsize <- c(2.5,2.5)
  
# Fold change is not a measure of change, it includes how much you have to start with
# FC = Y/X

# Define Conversion table between fc<1 and fc>=1
y_rats <- c(2, 3, 4, 5, 6)
df_fc = data.frame(nfc = 1/rev(y_rats), pfc = y_rats)

# Test dataframe for fc conversion functions
fc_test <- data.frame(ind = seq(1,2*nrow(df_fc)+1), 
                      fcu0 = seq(1,2*nrow(df_fc)+1) - length(y_rats) -1,
                      fc = c(df_fc$nfc, 1, df_fc$pfc))
fc_test$log2fc <- log2(fc_test$fc)
fc_test$mfc <- fc_to_mfc(fc_test$fc)
fc_test$con_mfc <- contract1(fc_test$mfc)

##  Linear Plot
g0 <- ggplot( data = fc_test, aes(y = fc, x = fcu0)) + 
  geom_point(size = 1) + 
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(breaks = fcu0,
                     labels = fcu0) +
  ylab("FC") + xlab("FC from No Change") +
  theme_minimal(base_size = 8) + theme(panel.grid.minor = element_blank())
g0
save_plot(paste(fig_path, '/', "A_fc_linear.jpg", sep = ""),
          g0, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


## Log plot of fold change
g1 <- ggplot( data = fc_test, aes(y = fc, x = fcu0)) + 
  geom_point(size = 1) + 
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 0) +
  scale_y_continuous(trans='log2',
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  scale_x_continuous(breaks = fcu0, labels = fcu0) +
  ylab("FC") + xlab("FC from No Change") +
  theme_minimal(base_size = 8) + theme(panel.grid.minor = element_blank())
g1
save_plot(paste(fig_path, '/', "B_fc_log2.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])

# log of fold change
g2 <- ggplot( data = fc_test, aes(y = log2fc, x = fcu0)) + 
  geom_point(size = 1) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(breaks = fcu0, labels = fcu0) +
  ylab("Log2 FC") + xlab("FC from No Change") +
  theme_minimal(base_size = 8) + theme(panel.grid.minor = element_blank())
g2
save_plot(paste(fig_path, '/', "B_log2_fc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])


# Mirrored Contracted Fold Change
g3 <- ggplot( data = fc_test, aes(y = con_mfc, x = fcu0)) + 
  geom_point() + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(breaks = fcu0, labels = fcu0) +
  ylab("Con-FC") + xlab("FC from No Change") +
  theme_minimal(base_size = 8) + theme(panel.grid.minor = element_blank())
g3
save_plot(paste(fig_path, '/', "C_con_mfc.jpg", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])


# Mirror fold change plot
g4 <- ggplot( data = fc_test, aes(y = con_mfc, x = fcu0)) + 
  geom_point() + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_y_continuous(breaks = fc_test$con_mfc, labels = fc_test$con_mfc) +
  scale_x_continuous(breaks = fcu0,    labels = fcu0) +
  ylab("MAD-FC") + xlab("FC from No Change") +
  theme_minimal(base_size = 8) + theme(panel.grid.minor = element_blank())
g4
  g4 <- gg_revaxis_mfc(g4,'y', num_format = "fraction")
g4
save_plot(paste(fig_path, '/', "C_mad-fc.jpg", sep = ""),
          g4, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])