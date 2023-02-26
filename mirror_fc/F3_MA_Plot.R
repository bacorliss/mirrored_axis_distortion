


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
install.packages(
  pkgs = "DESeqAnalysis",
  repos = c(
    "https://r.acidgenomics.com",
    BiocManager::repositories()
  ),
  dependencies = TRUE
)


# Linear visualization
base_dir = "mirror_fc"
fig_num = "3" 
fig_path = paste(base_dir,"/figure/F",fig_num, sep="")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
ggsize <- c(2.25,2.25)

library(ggplot2)
library(DESeqAnalysis)
library(scales)

data(deseq)
res <- results(deseq)

res$log2FoldChange


df <- data.frame(baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, 
                 pvalue<-res$pvalue,padj = res$padj)

ggplot(data=df,aes(x=baseMean, y=log2FoldChange, color = pvalue <0.05)) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 0.8, size = 1) +
  scale_x_continuous(trans='log10',labels = label_comma()) +
  theme_classic(base_size = 8)


## Get genes from DESeqDataSet.
dds <- as(deseq, "DESeqDataSet")
genes <- head(rownames(dds))
print(genes)
#> [1] "gene1" "gene2" "gene3" "gene4" "gene5" "gene6"

## DESeqAnalysis ====
plotMA(dds, i = 1L)