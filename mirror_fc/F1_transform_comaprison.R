

library(ggplot2)
library(cowplot)

source("R/mirrored_axis_distortion.R")

# Linear visualization
base_dir = "mirror_fc"
fig_num = "1" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)


# Fold change is not a measure of change, it includes how much you have to start with
# FC = Y/X

# Define Conversion table between fc<1 and fc>=1
y_rats <- c(2, 3, 4, 5, 6, 7)
df_fc = data.frame(nfc = 1/rev(y_rats), pfc = y_rats)

# Test dataframe for fc conversion functions
fc_test <- data.frame(ind = seq(1,2*nrow(df_fc)+1), fc = c(df_fc$nfc, 1, df_fc$pfc))
fc_test$log2fc <- log2(fc_test$fc)
fc_test$mfc <- fc_to_mfc(fc_test$fc)


ggplot( data = fc_test, aes(y = log2fc, x = fc)) + 
  geom_point() + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1) +
  # scale_y_continuous(trans = "log2") +
  ylab("Log2 FC") + xlab("FC") +
  theme_minimal()



ggplot( data = fc_test, aes(y = fc, x = fc)) + 
  geom_point() + 
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  ylab("FC") + xlab("FC") +
  theme_minimal()



ggplot( data = fc_test, aes(y = mfc, x = fc)) + 
  geom_point() + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1) +
  ylab("MFC") + xlab("FC") +
  theme_minimal()




# Define Conversion table between fc<1 and fc>=1
y_rats <- c(10.1/10, 11/10, 12.5/10, 15/10, 2, 3, 4, 5)
df_fc = data.frame(nfc = 1/y_rats, pfc = y_rats)

# Test dataframe for fc conversion functions
fc_test <- 
  data.frame(ind = c(1, 1+seq(1,nrow(df_fc)), 1+seq(1,nrow(df_fc))), Direction = c("NC", rep("-FC",nrow(df_fc)),
  rep("+FC",nrow(df_fc))), fc = c( 1, df_fc$nfc, df_fc$pfc))
fc_test$mfc <- mirror_fc(fc_test$fc)
fc_test$contract_mfc <- contract1(fc_test$mfc)


manual_colors = c("#d95f02", "#1b9e77", "#7570b3")




# fc_test$mfc_2_fc <- mirror_fc(fc_test$fc_2_mfc, forward = FALSE)
ggsize = c(2,2)

g1 <- ggplot(data = fc_test, aes(y = ind, x = fc, color = Direction)) + 
  geom_vline(xintercept = 1) + scale_y_reverse() +
  geom_point(size=1) + theme_minimal(base_size = 8) + ylab("Row") + xlab("FC") + 
  theme(legend.position="none", panel.grid.minor.x = element_blank()) +
  scale_y_continuous(labels = as.character(fc_test$ind), breaks = fc_test$ind) +
  scale_color_manual(values=manual_colors)
g1
save_plot(paste(fig_path, '/', "A_fc.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])


g2 <- ggplot(data = fc_test, aes(y = ind, x = mfc, color = Direction)) + 
  annotate("rect", ymax = Inf, ymin = -Inf, xmin = -1, xmax = 1, alpha = 0.25) +
  geom_point(size=1) + theme_minimal(base_size = 8) + ylab("Row") + xlab("MFC") + 
  theme(legend.position="none", panel.grid.minor.x = element_blank()) + scale_y_reverse() +
  scale_y_continuous(labels = as.character(fc_test$ind), breaks = fc_test$ind) +
  scale_color_manual(values=manual_colors) 
g2
save_plot(paste(fig_path, '/', "B_mfc.jpg", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])


g3 <- ggplot(data = fc_test, aes(y = ind, x = contract_mfc, color = Direction)) + 
  geom_vline(xintercept = 0) +
  geom_point(size=1) + theme_minimal(base_size = 8) + ylab("Row") + xlab("Con-MFC") + 
  theme(legend.position="none", panel.grid.minor.x = element_blank()) + scale_y_reverse() +
  scale_y_continuous(labels = as.character(fc_test$ind), breaks = fc_test$ind)  +
  scale_color_manual(values=manual_colors)
g3
save_plot(paste(fig_path, '/', "C_con-mfc.jpg", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])

g4 <- gg_revaxis_mfc(g3,'x',num_format = "decimal") + xlab("MAD-FC") 
g4
save_plot(paste(fig_path, '/', "D_Ax-Rev_con-mfc_decimal.jpg", sep = ""),
          g4, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])

g5 <- gg_revaxis_mfc(g3,'x',num_format = "power") + xlab("MAD-FC") 
g5
save_plot(paste(fig_path, '/', "D_Ax-Rev_con-mfc_power.jpg", sep = ""),
          g5, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])


g6 <- gg_revaxis_mfc(g3,'x',num_format = "fraction") + xlab("MAD-MFC") 
g6
save_plot(paste(fig_path, '/', "D_Ax-Rev_con-mfc_fraction.jpg", sep = ""),
          g6, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])


# scale_x_continuous(labels = new_xlabs, breaks = x_breaks)










