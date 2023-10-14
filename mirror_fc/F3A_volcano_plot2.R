

# Dataset: Airway package in bioconductor
# Citation: https://pubmed.ncbi.nlm.nih.gov/24926665/

# This code is appropriated from the introduction vignette for DESeq2
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


# Install required Base Packages
base_packages <- c("magrittr", "ggplot2", "tidyverse", "cowplot","BiocManager")
install.packages(setdiff(base_packages, rownames(installed.packages())))  
# Install required Bioconductor Packages
biocm_packages <- c("DESeq2", "org.Hs.eg.db","airway","EnhancedVolcano")
bioc_installs <- setdiff(biocm_packages, rownames(installed.packages()))
if (length(bioc_installs)) {BiocManager::install(bioc_installs) }

# Load base packages
lapply(base_packages, library, character.only = TRUE)
# Load Bioconductor packages packages
lapply(biocm_packages, library, character.only = TRUE)
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

# Calculate fold change
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds)
res <- results(dds,
               contrast = c('dex','trt','untrt'))
res <- lfcShrink(dds,
                 contrast = c('dex','trt','untrt'), res=res, type = 'normal')


df_res <- as.data.frame(res)
df_res$Change <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_res$Change[df_res$log2FoldChange > 1 & df_res$pvalue < 0.1] <- "Up"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_res$Change[df_res$log2FoldChange < -1 & df_res$pvalue < 0.1] <- "Down"
df_res$Change <- factor(df_res$Change, ordered = TRUE, levels = c("Up", "NS", "Down"))
df_res$FoldChange <- 2^df_res$log2FoldChange
df_res$mcFoldChange <- contract1(fc_to_mfc(df_res$FoldChange))



# Log2 Plot
g1 <- ggplot(data=df_res, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(col=Change), size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic(base_size = 8) + 
  ylab(expression(-log[10]~(`p-value`))) +
  xlab(expression(log[2]~(FC))) +
  theme(legend.position = "none")     
g1
save_plot(paste(fig_path, '/', "A_volcano_log2.png", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# Fold Change Linear Plot
g2 <- ggplot(data=df_res, aes(x=FoldChange, y=-log10(pvalue))) +
  geom_vline(xintercept = 1, color = "gray") +
  geom_point(aes(col=Change), size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic(base_size = 8)   +
  ylab(expression(-log[10]~(`p-value`))) +
  xlab("FC") +
  theme(legend.position = "none")     
g2
save_plot(paste(fig_path, '/', "B_volcano_FC.png", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  


# MAD Fold Change Linear Plot
g3 <- ggplot(data=df_res, aes(x=mcFoldChange, y=-log10(pvalue))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(col=Change), size = 0.5, alpha = 0.6, shape = 16) +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8) +
  ylab(expression(-log[10]~(`p-value`))) +
  xlab("FC") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(c(-9, - 4, 0, seq(4,30,5))), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))
g3
g3 <- gg_revaxis_mfc(g3,'x', num_format = "fraction")
g3
save_plot(paste(fig_path, '/', "C_volcano_MFC.png", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])  

