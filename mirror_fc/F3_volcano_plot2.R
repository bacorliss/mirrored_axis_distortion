

# Citation: https://pubmed.ncbi.nlm.nih.gov/24926665/

# https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}

if (!requireNamespace('EnhancedVolcano', quietly = TRUE)) {
  BiocManager::install('EnhancedVolcano')
}; library(EnhancedVolcano)

library(airway)
library(magrittr)
library(org.Hs.eg.db)
library('DESeq2')


source("R/mirrored_axis_distortion.R")




data('airway')
airway$dex %<>% relevel('untrt')

ens <- rownames(airway)


symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(airway), names(symbols))]
rownames(airway) <- symbols
keep <- !is.na(rownames(airway))
airway <- airway[keep,]



dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds,
               contrast = c('dex','trt','untrt'))
res <- lfcShrink(dds,
                 contrast = c('dex','trt','untrt'), res=res, type = 'normal')
res$FoldChange <- 2^res$log2FoldChange



library(plotly)

g<-EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

ggplotly(g)

df_gene <- data.frame( name = rownames(res), log2FoldChange = res$log2FoldChange, 
                       FoldChange = 2^res$log2FoldChange, pvalue = res$pvalue)


df_gene$diffexpressed <- ""
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df_gene$diffexpressed[df_gene$log2FoldChange > 0.5 & df_gene$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df_gene$diffexpressed[df_gene$log2FoldChange < -0.5 & df_gene$pvalue < 0.05] <- "DOWN"





df_gene$mcFoldChange <- contract1(fc_to_mfc(df_gene$FoldChange))

# Log2 Plot
g1 <- ggplot(data=df_gene, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(col=diffexpressed)) +
  scale_color_manual(values=c("blue", "black", "red")) +
  theme_minimal() 
g1


# Fold Change Linear Plot
g2 <- ggplot(data=df_gene, aes(x=FoldChange, y=-log10(pvalue))) +
  geom_point(aes(col=diffexpressed)) +
  scale_color_manual(values=c("blue", "black", "red")) +
  theme_minimal()  
g2


# Fold Change Linear Plot
g3 <- ggplot(data=df_gene, aes(x=mcFoldChange, y=-log10(pvalue))) +
  geom_point(aes(col=diffexpressed)) +
  scale_color_manual(values=c("blue", "black", "red")) +
  theme_minimal() 
g3 <- gg_revaxis_mfc(g3,'x', num_format = "fraction")
g3

