

# Citation: https://pubmed.ncbi.nlm.nih.gov/24926665/

# https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}

if (!requireNamespace('EnhancedVolcano', quietly = TRUE)) {
  BiocManager::install('EnhancedVolcano')
}







library(airway)
library(magrittr)

data('airway')
airway$dex %<>% relevel('untrt')

ens <- rownames(airway)

library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(airway), names(symbols))]
rownames(airway) <- symbols
keep <- !is.na(rownames(airway))
airway <- airway[keep,]

library('DESeq2')

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds,
               contrast = c('dex','trt','untrt'))
res <- lfcShrink(dds,
                 contrast = c('dex','trt','untrt'), res=res, type = 'normal')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

