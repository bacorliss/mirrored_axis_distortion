# F3_FoldChange_Protein_expression

#https://pubmed.ncbi.nlm.nih.gov/25612154/

library(tidyverse)
source("R/mirrored_axis_distortion.R")

base_dir = "mirror_fc"
fig_num = "3" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)




df0 <- read.csv(paste0(base_dir, "/gene_expression_intervals.csv"))
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
ggplot(data = df,aes(x=label, y = log2_mid, color = color)) +
  geom_errorbar(aes(ymax = log2_hi, ymin = log2_lo), size = 1, width = .35,
                position=position_dodge(width=0.4)) + 
  geom_point(size = 2, position=position_dodge(width=0.4)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none") +
  scale_color_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c("#1b9e77", "#d95f02", "#7570b3")) +
  scale_shape_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c(18, 15, 17)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("") + ylab("Log2(FC) From Vehicle TX")


# FC Plot
ggplot(data = df,aes(x=label, y = mid, color = color, shape = color)) +
  geom_errorbar(aes(ymax = hi, ymin = lo), size = 1, width = .35,
                position=position_dodge(width=0.4)) + 
  geom_point(size = 2, position=position_dodge(width=0.4)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none") +
  scale_color_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c("#1b9e77", "#d95f02", "#7570b3")) +
  scale_shape_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c(18, 15, 17)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("") + ylab("Fold-change From Vehicle TX")



# MAD FC Plot
ggplot(data = df,aes(x=label, y = mfc_mid, color = color, shape = color)) +
  geom_errorbar(aes(ymax = mfc_hi, ymin = mfc_lo), size = 1, width = .35,
                position=position_dodge(width=0.4)) + 
  geom_point(size = 2, position=position_dodge(width=0.4)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none") +
  scale_color_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c("#1b9e77", "#d95f02", "#7570b3")) +
  scale_shape_manual(name = "Treatment", labels = c("MEKi", "PI3Ki", "MEKI:PI3Ki"), values = c(18, 15, 17)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("") + ylab("Fold-change From Vehicle TX")
