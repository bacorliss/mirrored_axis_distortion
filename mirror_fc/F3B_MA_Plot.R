

# Code appropriated from:
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#the-deseqdataset


# Install required Base Packages
base_packages <- c("magrittr", "ggplot2", "tidyverse", "cowplot","BiocManager")
install.packages(setdiff(base_packages, rownames(installed.packages())))  
# Install required Bioconductor Packages
biocm_packages <- c("DESeq2", "org.Hs.eg.db","airway")
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
ggsize <- c(2.25,2.15)


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
# Perform differential gene expression analysis with Walde significance test without a beta prior
dds <- DESeq(dds, betaPrior=FALSE)
# Extract results with treated versus untreated cells 
res <- results(dds,
               contrast = c('dex','trt','untrt'))
# Perform shrinkage of dispersion estimates to refine log2 fold changes estimates
res <- lfcShrink(dds,
                 contrast = c('dex','trt','untrt'), res=res, type = 'normal')
# Fold changes were calcualted from 
res$FoldChange <- 2^res$log2FoldChange



# Convert DESeq data to dataframe for plotting with ggplot
df_gene <- data.frame( name = rownames(res), baseMean = res$baseMean, 
                       log2FoldChange = res$log2FoldChange, 
                       FoldChange = 2^res$log2FoldChange, pvalue = res$pvalue)
df_gene$diffexpressed <- "None"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_gene$diffexpressed[df_gene$log2FoldChange > 1 & df_gene$pvalue < 0.05] <- "Up"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_gene$diffexpressed[df_gene$log2FoldChange < -1 & df_gene$pvalue < 0.05] <- "Down"

df_gene$diffexpressed <- factor(df_gene$diffexpressed, ordered = TRUE, levels = c("Up", "None", "Down"))
df_gene$mcFoldChange <- contract1(fc_to_mfc(df_gene$FoldChange))





# Log2 FC Plot
g1 <- ggplot(data = df_gene, aes(x = baseMean, y = log2FoldChange)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_point(aes(col=diffexpressed), size = 0.4, alpha = 0.6, shape = 16) +
  scale_x_log10() +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  theme_classic(base_size = 8) + 
  ylab(expression(log[2]~(FC))) +
  xlab("Normalized Mean Count") +
  theme(legend.position = "none")     
g1
save_plot(paste(fig_path, '/', "A_ma_log2_plot.png", sep = ""),
          g1, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])


# Linear FC Plot
g2 <- ggplot(data = df_gene, aes(x = baseMean, y = FoldChange)) +
  geom_hline(yintercept = 1, color = "grey") +
  geom_point(aes(col=diffexpressed), size = 0.4, alpha = 0.6, shape = 16) +
  scale_x_log10() +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  theme_classic(base_size = 8) + 
  ylab("FC") +
  xlab("Normalized Mean Count") +
  theme(legend.position = "none")     
g2
save_plot(paste(fig_path, '/', "B_ma_fc_plot.png", sep = ""),
          g2, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])


# MAD-FC Plot
g3 <- ggplot(data = df_gene, aes(x = baseMean, y = mcFoldChange)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_point(aes(col=diffexpressed), size = 0.4, alpha = 0.6, shape = 16) +
  scale_x_log10() +
  scale_color_manual(values=c("blue", "black", "red"), name = "Change") +
  scale_y_continuous(breaks = c(-10, 0, 10, 20)) +
  theme_classic(base_size = 8) + 
  ylab("MAD-FC") +
  xlab("Normalized Mean Count") +
  theme(legend.position = "none")     
g3
g3 <- gg_revaxis_mfc(g3,'y', num_format = "fraction")
g3

save_plot(paste(fig_path, '/', "C_ma_mad-fc_plot.png", sep = ""),
          g3, dpi = 600, base_height = ggsize[1], 
          base_width = ggsize[2])


