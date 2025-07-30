############### packages
library(dplyr)
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

source("GEM_BXD_CECUM_functions_zy.R")



############### load data 
metadata <- read.table("S1_metadata_summary_zy.txt", header=T, row.names= 1, sep="\t") 
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

MG_pre <- BXD.Data.Pre(df = MGabun, datatype = "MG", taxlevel = "species")
BXD.Anova.plot(MG_pre, metadata)

MT_pre <- BXD.Data.Pre(df = MTabun, datatype = "MT", taxlevel = "species")
BXD.Anova.plot(MT_pre, metadata)



############### Fig1 E&F 
degDiet <- BXD.DEG.plot(mRNAraw, metadata, label_n = 20, group = "Diet")
degAge <- BXD.DEG.plot(mRNAraw, metadata, label_n = 20, group = "Age")


############### Fig1 H 


############### Fig2 A 


############### Fig2 B&C 


############### Fig2 D &  Fig3 A



############### Fig1 E-G 
### https://github.com/biobakery/maaslin3


############### Fig3 B 
### https://github.com/biobakery/humann

############### Fig3 C&D 






