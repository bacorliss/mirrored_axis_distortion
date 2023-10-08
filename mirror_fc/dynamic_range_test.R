

library(ggplot2)
library(cowplot)
library(patchwork)


# Linear visualization
base_dir = "mirror_fc"
fig_num = "1" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(3,6)




colpal = c("black", "#FA6454", "black", "#95B8D7", "black", "#CB86F6", "black", "#32DE6B","black", "#FF9900")



# Dataset 1: 3 orders of magnitude
df <- data.frame(x=c(1, 2,2^2, 2^3, 2^4, 2^5,2^6,2^7,2^8,2^9),
                 y=c(rep(c(-1,-.75),5)),
                 yend = c(rep(c(0.75,1),5)))
df$Color <- factor(rep(1:length(colpal), (nrow(df) %/% length(colpal))+1) [1:nrow(df)])


mycolors <- rep(colpal, (nrow(df) %/% length(colpal))+1) [1:nrow(df)]

p1 <- ggplot(df, aes(x=x, y=y, color = Color)) +
  geom_segment(aes(x=x,xend=x,y=y,yend=yend), size =1.2) +
  scale_color_manual(values=mycolors)+
  labs(x="",y="") +
  theme_classic() + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.position="none")

# save_plot(paste(fig_path, '/', "dynamic_range.jpg", sep = ""),
#           g1, dpi = 600, base_height = ggsize[1], 
#           base_width = ggsize[2])


p2 <- ggplot(df, aes(x=x, y=y, color = Color)) +
  geom_segment(aes(x=x,xend=x,y=y,yend=yend), size =1.2) +
  scale_x_continuous(trans='log2') +
  scale_color_manual(values=mycolors)+
  labs(x="",y="") +
  theme_classic() + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.position="none")


g1 <-p1 + plot_spacer() + p2 + plot_layout(ncol = 1, heights=c(2, -.5 , 2))
g1


# g1 <- plot_grid(p1, p2, ncol = 1, align = "v")

save_plot(paste(fig_path, '/', "dynamic_range.jpg", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])




