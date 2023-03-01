

# Citation: https://pubmed.ncbi.nlm.nih.gov/24926665/

# https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}

if (!requireNamespace('EnhancedVolcano', quietly = TRUE)) {
  BiocManager::install('EnhancedVolcano')
}; library(EnhancedVolcano)

library('airway')
library('magrittr')
library('org.Hs.eg.db')
library('DESeq2')
library("cowplot")
source("R/mirrored_axis_distortion.R")


# Linear visualization
base_dir = "mirror_fc"
fig_num = "3" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2.25,2.25)


data('airway')
airway$dex %<>% relevel('untrt')
ens <- rownames(airway)

# Map symbols
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(airway), names(symbols))]
rownames(airway) <- symbols
keep <- !is.na(rownames(airway))
airway <- airway[keep,]

# Calcualte fold change
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds,
               contrast = c('dex','trt','untrt'))
res <- lfcShrink(dds,
                 contrast = c('dex','trt','untrt'), res=res, type = 'normal')
res$FoldChange <- 2^res$log2FoldChange







df_gene <- data.frame( name = rownames(res), log2FoldChange = res$log2FoldChange, 
                       FoldChange = 2^res$log2FoldChange, pvalue = res$pvalue)
df_gene$diffexpressed <- "None"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_gene$diffexpressed[df_gene$log2FoldChange > 1 & df_gene$pvalue < 0.05] <- "Up"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_gene$diffexpressed[df_gene$log2FoldChange < -1 & df_gene$pvalue < 0.05] <- "Down"

df_gene$diffexpressed <- factor(df_gene$diffexpressed, ordered = TRUE, levels = c("Up", "None", "Down"))
df_gene$mcFoldChange <- contract1(fc_to_mfc(df_gene$FoldChange))



# Log2 Plot
g1 <- ggplot(data=df_gene, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(aes(col=diffexpressed), size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8) + 
  ylab(expression(-log[10]~(p-value))) +
  xlab(expression(log[2]~(FC))) +
  theme(legend.position = "none")     
g1
save_plot(paste(fig_path, '/', "A_volcano_log2.png", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# Fold Change Linear Plot
g2 <- ggplot(data=df_gene, aes(x=FoldChange, y=-log10(pvalue))) +
  geom_vline(xintercept = 1, color = "grey") +
  geom_point(aes(col=diffexpressed), size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8)   +
  ylab(expression(-log[10]~(p-value))) +
  xlab("FC") +
  theme(legend.position = "none")     
g2
save_plot(paste(fig_path, '/', "B_volcano_FC.png", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# MAD Fold Change Linear Plot
g3 <- ggplot(data=df_gene, aes(x=mcFoldChange, y=-log10(pvalue))) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(aes(col=diffexpressed), size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8) +
  ylab(expression(-log[10]~(p-value))) +
  xlab("MAD-FC") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(c(-9, - 4, 0, seq(4,30,5))), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))
g3
g3 <- gg_revaxis_mfc(g3,'x', num_format = "fraction")
g3
save_plot(paste(fig_path, '/', "C_volcano_MFC.png", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

