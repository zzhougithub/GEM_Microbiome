############### packages - These are all the libraries 
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
library(readxl)
library(maaslin3)
library(knitr)

source("GEM_BXD_CECUM_functions_zy.R")



############### load data 
metadata <- read_excel("S1_metadata_cecalsamples.xlsx", sheet = 2, na = character()) %>%
  column_to_rownames(var = names(.)[1])
mRNAraw <- read.table("S2_mRNA_cecum_raw_counts_zy.txt", header=T, sep="\t") 
mRNAlog2 <- read.table("S3_mRNA_cecum_TPMLog2_zy.txt", header=T, sep="\t") 
MGabun <- read.table("S4_metaDNA_metaphlan4_abundance_zy.txt", header=T, na.strings="", sep="\t", stringsAsFactors=FALSE)
MTabun <- read.table("S5_metaRNA_metaphlan4_abundance_zy.txt", header=T, na.strings="", sep="\t", stringsAsFactors=FALSE)
MGpath <- read_tsv("S6_metaDNA_humann3.8_pathabundance_relab_zy.txt", col_types = cols(), na = "")
MTpath <- read_tsv("S7_metaRNA_humann3.8_pathabundance_relab_zy.txt", col_types = cols(), na = "")

## For the "metadata" data frame, take care that some metadata analyses require removing replicate samples (all those with _r or _n suffixes)!! These are technical replicates of the exact same mouse.
  
  
############### Fig1 A&B ; 
# Hand-drawn figures

############### Fig1 C ; FigS1 D 
mRNA_pre <- BXD.Data.Pre(df = mRNAlog2, datatype = "mRNA")
BXD.Anova.plot(mRNA_pre, metadata)
BXD.Fstat.plot(mRNA_pre, metadata)

MG_pre <- BXD.Data.Pre(df = MGabun, datatype = "abundance", taxlevel = "species")
BXD.Anova.plot(MG_pre, metadata)

MT_pre <- BXD.Data.Pre(df = MTabun, datatype = "abundance", taxlevel = "species")
BXD.Anova.plot(MT_pre, metadata)

############### Fig1 D ; 

# This code is in the separate file 2_VariancesExplained_Code.r

############### Fig S1 F

# This code is in the separate file 2_MGMTCorrelations_FigS1F.R


############### Fig1 E&F 
degDiet <- BXD.DEG.plot(mRNAraw, metadata, label_n = 20, group = "Diet")
degAge <- BXD.DEG.plot(mRNAraw, metadata, label_n = 20, group = "Age")


############### Fig1 G, Left (within liver data)
# Uses the below code, but on the liver data from Supplemental Data 1 from our previous publication on the same mice but looking at liver gene expression data
# https://www.sciencedirect.com/science/article/pii/S2405471221003446?via%3Dihub#mmc8

# CD and HF variables use the following row:
# mRNA_43016_Cyp2c55

## !! NOTE THAT THIS DATA IS FROM THE PREVIOUS STUDY AND IS NOT INCLUDED AS SUPPLEMENTAL IN THIS PAPER!!
LiverDat_header=read.table("aData_S1_AllOmicsandPhenotypeData.csv", nrows=8, row.names=1, header=T, sep=",")
LiverDat_header=LiverDat_header[,3:ncol(LiverDat_header)]

liverDatlog2=read.table("aData_S1_AllOmicsandPhenotypeData.csv", row.names=1, header=T,  na.strings="", sep=",", skip=8, stringsAsFactors=FALSE)
liverDatlog2=liverDatlog2[1:nrow(liverDatlog2),3:ncol(liverDatlog2)]

livermRNAlog2Diets=as.character(LiverDat_header['Diet',])

CD=as.numeric(liverDatlog2['mRNA_43016_Cyp2c55',])
HF=as.numeric(liverDatlog2['mRNA_43016_Cyp2c55',])

CD <- CD[livermRNAlog2Diets == "CD"]
HF <- HF[livermRNAlog2Diets == "HF"]

legendposition="topleft"

groupA=rm.outlier(CD)
groupA=rm.outlier(groupA)

groupB=rm.outlier(HF)
stripchart(groupA, at=0.7, pch=21, cex=1.2, bg="lightblue", method="jitter", xlim=c(0.6, 1.3), las=1, ylim=c(0,2.8), vertical=T)
stripchart(groupB, at=1.0, pch=21, cex=1.2, bg="orange", method="jitter", add=T, vertical=T)
	segments(0.65, quantile(groupA, na.rm=T)[2], 0.75, quantile(groupA, na.rm=T)[2], col="red", lwd=3)	# first quantile bar for first cohort
	segments(0.6, median(groupA, na.rm=T), 0.8, median(groupA, na.rm=T), col="red", lwd=3)	# median bar for first cohort
	segments(0.65, quantile(groupA, na.rm=T)[4], 0.75, quantile(groupA, na.rm=T)[4], col="red", lwd=3)	# first quantile bar for first cohort
	segments(0.7,quantile(groupA, na.rm=T)[2],0.7,quantile(groupA, na.rm=T)[4], col="red", lwd=3)
	segments(.95, quantile(groupB, na.rm=T)[2], 1.05, quantile(groupB, na.rm=T)[2], col="red", lwd=3)	# first quantile bar for second cohort	
	segments(.9, median(groupB, na.rm=T), 1.1, median(groupB, na.rm=T), col="red", lwd=3)	# median bar for second cohort, etc
	segments(.95, quantile(groupB, na.rm=T)[4], 1.05, quantile(groupB, na.rm=T)[4], col="red", lwd=3)	# first quantile bar for second cohort
	segments(1.0,quantile(groupB, na.rm=T)[2],1.0,quantile(groupB, na.rm=T)[4], col="red", lwd=3)
# legend	
pr = signif(t.test(groupA,groupB, paired=FALSE)$p.value,digits=2)			# This outputs the p-value of the correlation
pletter=substitute(paste("   ", italic("p"), " = ", pr), list(pr=signif(pr[1], digits=2)))	# Makes the r-value with italicized r

legend(legendposition, bg='white', bty="n", cex=1, text.col="red", legend=c(as.expression(pletter)))



############### Fig1 G, Right

mRNAlog2Diets=c("HF","CD","HF","CD","HF","HF","CD","HF","CD","CD","CD","HF","HF","HF","CD","HF","CD","CD","CD","HF","HF","CD","CD","CD","HF","HF","CD","HF","HF","HF","HF","HF","HF","HF","HF","HF","CD","CD","CD","CD","CD","HF","HF","HF","CD","HF","CD","CD","HF","HF","CD","CD","CD","HF","CD","HF","HF","CD","HF","HF","HF","CD","CD","CD","CD","CD","HF","CD","CD","HF","HF","CD","HF","HF","CD","CD","HF","CD","HF","HF")
mRNAlog2 <- read.table("S3_mRNA_cecum_TPMLog2_zy.txt", row.names=1, header=T, sep="\t") 

CD=as.numeric(mRNAlog2['Cyp2c55',])
HF=as.numeric(mRNAlog2['Cyp2c55',])

CD <- CD[mRNAlog2Diets == "CD"]
HF <- HF[mRNAlog2Diets == "HF"]

legendposition="topleft"

groupA=rm.outlier(CD)
groupA=rm.outlier(groupA)

groupB=rm.outlier(HF)
stripchart(groupA, at=0.7, pch=21, cex=1.2, bg="lightblue", method="jitter", xlim=c(0.6, 1.3), las=1, ylim=c(10.8,14), vertical=T)
stripchart(groupB, at=1.0, pch=21, cex=1.2, bg="orange", method="jitter", add=T, vertical=T)
	segments(0.65, quantile(groupA, na.rm=T)[2], 0.75, quantile(groupA, na.rm=T)[2], col="red", lwd=3)	# first quantile bar for first cohort
	segments(0.6, median(groupA, na.rm=T), 0.8, median(groupA, na.rm=T), col="red", lwd=3)	# median bar for first cohort
	segments(0.65, quantile(groupA, na.rm=T)[4], 0.75, quantile(groupA, na.rm=T)[4], col="red", lwd=3)	# first quantile bar for first cohort
	segments(0.7,quantile(groupA, na.rm=T)[2],0.7,quantile(groupA, na.rm=T)[4], col="red", lwd=3)
	segments(.95, quantile(groupB, na.rm=T)[2], 1.05, quantile(groupB, na.rm=T)[2], col="red", lwd=3)	# first quantile bar for second cohort	
	segments(.9, median(groupB, na.rm=T), 1.1, median(groupB, na.rm=T), col="red", lwd=3)	# median bar for second cohort, etc
	segments(.95, quantile(groupB, na.rm=T)[4], 1.05, quantile(groupB, na.rm=T)[4], col="red", lwd=3)	# first quantile bar for second cohort
	segments(1.0,quantile(groupB, na.rm=T)[2],1.0,quantile(groupB, na.rm=T)[4], col="red", lwd=3)
# legend	
pr = signif(t.test(groupA,groupB, paired=FALSE)$p.value,digits=2)			# This outputs the p-value of the correlation
pletter=substitute(paste("   ", italic("p"), " = ", pr), list(pr=signif(pr[1], digits=2)))	# Makes the r-value with italicized r

legend(legendposition, bg='white', bty="n", cex=1, text.col="red", legend=c(as.expression(pletter)))


############### Fig1 H  ; FigS2 A
BXD.ORA.plot(degAgeHF, 15)
BXD.ORA.plot(degDiet, 15)


############### Fig2 B&C 
MG_phyl <- BXD.Data.Pre(df = MGabun, datatype = "abundance", taxlevel = "phylum")
BXD.Ratio.plot(MG_phyl, metadata, "p__Firmicutes", "p__Bacteroidetes")


############### Fig2 D
MG_genu <- BXD.Data.Pre(df = MGabun, datatype = "abundance", taxlevel = "genus")
MT_genu <- BXD.Data.Pre(df = MTabun, datatype = "abundance", taxlevel = "genus")
BXD.Stacked.plot(MG_genu, MT_genu, metadata, "taxonomy", "MG_", "MT_")

############### Fig2 E
### https://github.com/biobakery/maaslin3
metadata$Diet <- factor(metadata$Diet, levels = c("CD", "HF"))

age_idx <- rownames(metadata)[metadata$Age < 1000]
age_idx <- intersect(age_idx, colnames(MG_pre))
feature_table <- MG_pre[, age_idx]

output_dir <- "/your/output/folder"

set.seed(1)
fit_out <- maaslin3(input_data = feature_table,
                    input_metadata = metadata,
                    output = output_dir,
                    min_abundance = 0.0001,
                    min_prevalence = 0.1,
                    formula = '~Diet + (1|Age) + (1|Body_weight) + (1|Strain) + (1|Sex)',
                    reference = c("Diet,CD"),
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.25,
                    correction = "BH",
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    coef_plot_vars = c('Diet HF'), 
                    heatmap_vars = c('Diet HF'), 
                    max_pngs = 250,
                    cores = 4)
                    
                    
                    
###### Fig3 A
### All samples were at least 50% other, so the y-axis scale was modified to have a broken axis above 50% so the lower part of the graph was more visible.

MGpath_pre <- BXD.Data.Pre(df = MGpath, datatype = "pathabun")
MTpath_pre <- BXD.Data.Pre(df = MTpath, datatype = "pathabun")
BXD.Stacked.plot(MGpath_pre, MTpath_pre, metadata, "Pathway", "MG_", "MT_")













