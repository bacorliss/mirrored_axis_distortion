


# # Preprocessing code appropriated from:
# # https://hbctraining.github.io/Intro-to-R-with-DGE/lessons/B1_DGE_visualizing_results.html
# 
# ### Set thresholds
# padj.cutoff <- 0.05
# lfc.cutoff <- 0.58
# 
# 
# threshold <- res_tableOE$padj < padj.cutoff & abs(res_tableOE$log2FoldChange) > lfc.cutoff
# 
# 
# res_tableOE$threshold <- threshold      
# 
# sigOE <- data.frame(subset(res_tableOE, threshold==TRUE))


# https://github.com/mousepixels/sanbomics_scripts/blob/main/tutorial_complex_Heatmap.Rmds




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
fig_num = "7" 
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
dds0 <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds0, betaPrior=FALSE)
# res <- results(dds,
#                contrast = c('dex','trt','untrt'))
# res <- lfcShrink(dds,
#                  contrast = c('dex','trt','untrt'), res=res, type = 'normal')
# res$FoldChange <- 2^res$log2FoldChange


#as.data.frame(colData(airway))
# SampleName    cell   dex albut        Run avgLength Experiment    Sample    BioSample
# SRR1039508 GSM1275862  N61311 untrt untrt SRR1039508       126  SRX384345 SRS508568 SAMN02422669
# SRR1039509 GSM1275863  N61311   trt untrt SRR1039509       126  SRX384346 SRS508567 SAMN02422675
# SRR1039512 GSM1275866 N052611 untrt untrt SRR1039512       126  SRX384349 SRS508571 SAMN02422678
# SRR1039513 GSM1275867 N052611   trt untrt SRR1039513        87  SRX384350 SRS508572 SAMN02422670
# SRR1039516 GSM1275870 N080611 untrt untrt SRR1039516       120  SRX384353 SRS508575 SAMN02422682
# SRR1039517 GSM1275871 N080611   trt untrt SRR1039517       126  SRX384354 SRS508576 SAMN02422673
# SRR1039520 GSM1275874 N061011 untrt untrt SRR1039520       101  SRX384357 SRS508579 SAMN02422683
# SRR1039521 GSM1275875 N061011   trt untrt SRR1039521        98  SRX384358 SRS508580 SAMN02422677


rlog_out <- rlog(dds0, blind=FALSE) #get normalized count data from dds object
# Get sample data from normalized counts
mx_raw<-assay(rlog_out)

mx_samples <- mx_raw[apply(mx_raw,1,function(x) all(x!=0)),]



mx_fc <- mx_samples[,seq(2,8,2)]/mx_samples[,seq(1,7,2)]
mx_fc[is.nan(mx_fc)]

mx_raw[apply(mx_raw,1,function(x) all(is.nan(x))),]



library(gplots)
heatmap.2(mx_samples, dendrogram="none")  