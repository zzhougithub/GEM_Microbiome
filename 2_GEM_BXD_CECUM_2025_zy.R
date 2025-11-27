############### packages
library(dplyr)
library(tidyr) 
library(tidyverse)
library(outliers)
library(gplots)
library(ggplot2)
library(vioplot)
library(reshape2)
library(DESeq2)
library(EnhancedVolcano)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(ggh4x)
library(ggtext)
library(ggalluvial)
library(ggpubr)
library(lemon)
library(ggsci)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(enrichplot)
library(DOSE)
library(glue)
library(scales)
library(GseaVis)

source("GEM_BXD_CECUM_functions_zy.R")



############### load data 
metadata <- read.table("S1.2_metadata_summary_zy.txt", header=T, row.names= 1, sep="\t") 
mRNAraw <- read.table("S2_mRNA_cecum_raw_counts_zy.txt", header=T, sep="\t") 
mRNAlog2 <- read.table("S3_mRNA_cecum_TPMLog2_zy.txt", header=T, sep="\t") 
MGabun <- read.table("S4_metaDNA_metaphlan4_abundance_zy.txt", header=T, na.strings="", sep="\t", stringsAsFactors=FALSE)
MTabun <- read.table("S5_metaRNA_metaphlan4_abundance_zy.txt", header=T, na.strings="", sep="\t", stringsAsFactors=FALSE)
MGpath <- read_tsv("S6_metaDNA_humann3.8_pathabundance_relab_zy.txt", col_types = cols(), na = "")
MTpath <- read_tsv("S7_metaRNA_humann3.8_pathabundance_relab_zy.txt", col_types = cols(), na = "")



############### Fig1 C&D ; FigS1 D 
mRNA_pre <- BXD.Data.Pre(df = mRNAlog2, datatype = "mRNA")
BXD.Anova.plot(mRNA_pre, metadata)
BXD.Fstat.plot(mRNA_pre, metadata)

MG_pre <- BXD.Data.Pre(df = MGabun, datatype = "abundance", taxlevel = "species")
BXD.Anova.plot(MG_pre, metadata)

MT_pre <- BXD.Data.Pre(df = MTabun, datatype = "abundance", taxlevel = "species")
BXD.Anova.plot(MT_pre, metadata)



############### Fig1 E&F 
degDiet <- BXD.DEG.plot(mRNAraw, metadata, label_n = 20, group = "Diet")
degAge <- BXD.DEG.plot(mRNAraw, metadata, label_n = 20, group = "Age")

#degDiet$signif_rank <- ifelse(degDiet$group == "NS", 1, 0)
#degDiet_sorted <- degDiet[order(degDiet$signif_rank, 
#                                -degDiet$log10p,
#                                decreasing = FALSE), ]
#degDiet_sorted$signif_rank <- NULL
#write.table(degDiet_sorted,
#            file = "degDiet_sorted.txt",
#            sep = "\t",
#            quote = FALSE,
#            col.names = NA)


############### Fig1 H  ; FigS2 A
BXD.ORA.plot(degAgeHF, 15)
BXD.ORA.plot(degDiet, 15)


############### Fig2 B&C 
MG_phyl <- BXD.Data.Pre(df = MGabun, datatype = "abundance", taxlevel = "phylum")
BXD.Ratio.plot(MG_phyl, metadata, "p__Firmicutes", "p__Bacteroidetes")


############### Fig2 D &  Fig3 A
MG_genu <- BXD.Data.Pre(df = MGabun, datatype = "abundance", taxlevel = "genus")
MT_genu <- BXD.Data.Pre(df = MTabun, datatype = "abundance", taxlevel = "genus")
BXD.Stacked.plot(MG_genu, MT_genu, metadata, "taxonomy", "MG_", "MT_")

MGpath_pre <- BXD.Data.Pre(df = MGpath, datatype = "pathabun")
MTpath_pre <- BXD.Data.Pre(df = MTpath, datatype = "pathabun")
BXD.Stacked.plot(MGpath_pre, MTpath_pre, metadata, "Pathway", "MG_", "MT_")


############### Fig1 E-G 
### https://github.com/biobakery/maaslin3


############### Fig3 B 
### https://github.com/biobakery/humann






