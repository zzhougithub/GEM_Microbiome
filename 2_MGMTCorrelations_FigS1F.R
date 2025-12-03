library(gplots)
library(outliers)
library(readxl)
setwd("/Users/evan.williams/Dropbox/OwnCloud/zy_Shared_Folder/Data/10_EvanAnalysis_Mar2024/")

################ FIGURE S1F, COMPARE CORRELATIONS FOR SAMPLES WITH MG AND MT DATA AND SAME EXACT TAXON, BY DIFFERENT RANK
### Note that in the figure itself we only show the genus-level correlation, but the below code calculates
### the correlation at any selected taxonomic rank - but since most analyses are at genus level that is what we use
### and also because genus-level taxonomic assignments are relatively canonical, but species-level taxonomic assignments are less definitive

data_in_DNA <- read.table("/Users/evan.williams/Dropbox/OwnCloud/zy_Shared_Folder/Data/10_EvanAnalysis_Mar2024/S4_metaDNA_metaphlan4_abundance_zy.txt",  header=T, row.names=1, sep="\t")
data_in_RNA <- read.table("/Users/evan.williams/Dropbox/OwnCloud/zy_Shared_Folder/Data/10_EvanAnalysis_Mar2024/S5_metaRNA_metaphlan4_abundance_zy.txt",  header=T, row.names=1, sep="\t")
data_in_pheno <- read_excel("/Users/evan.williams/Dropbox/OwnCloud/zy_Shared_Folder/Data/10_EvanAnalysis_Mar2024/S1_metadata_cecalsamples.xlsx",  header=T, row.names=1, sep=",", skip=11)

# This is just the columns (i.e. samples) that are measured at the phenotype and MG or MT level
OverlapDNA=data_in_DNA[,intersect(colnames(data_in_DNA), colnames(data_in_pheno))]
OverlapPhenoDNA=data_in_pheno[,intersect(colnames(data_in_pheno), colnames(data_in_DNA))]

OverlapRNA=data_in_RNA[,intersect(colnames(data_in_RNA), colnames(data_in_pheno))]
OverlapPhenoRNA=data_in_pheno[,intersect(colnames(data_in_pheno), colnames(data_in_RNA))]

OverlapDNA=OverlapDNA[-1,] # Removes the sum total "Bacteria"
OverlapRNA=OverlapRNA[-1,]

############### So now let's split the data by diet
header_in_pheno <- read.table("/Users/evan.williams/Dropbox/OwnCloud/zy_Shared_Folder/Data/10_EvanAnalysis_Mar2024/0_PHENOTYPES.csv",  header=T, row.names=1, sep=",", skip=6)

OverlapHeaderDNA=header_in_pheno[,intersect(colnames(header_in_pheno), colnames(data_in_DNA))]
OverlapHeaderRNA=header_in_pheno[,intersect(colnames(header_in_pheno), colnames(data_in_RNA))]

OverlapDNA_CD=matrix(NA, nrow=nrow(OverlapDNA), ncol=length(OverlapDNA))
for(i in 1:length(OverlapDNA)) {
	if(as.character(OverlapHeaderDNA["Diet",i])=="CD") {   OverlapDNA_CD[,i]=as.numeric(as.character(OverlapDNA[,i]))	}
} 
rownames(OverlapDNA_CD)=rownames(OverlapDNA)
colnames(OverlapDNA_CD)=colnames(OverlapDNA)

OverlapDNA_HF=matrix(NA, nrow=nrow(OverlapDNA), ncol=length(OverlapDNA))
for(i in 1:length(OverlapDNA)) {
	if(as.character(OverlapHeaderDNA["Diet",i])=="HF") {   OverlapDNA_HF[,i]=as.numeric(as.character(OverlapDNA[,i]))	}
} 
rownames(OverlapDNA_HF)=rownames(OverlapDNA)
colnames(OverlapDNA_HF)=colnames(OverlapDNA)

#####################
################ THE SAME COMPARISON FOR MT DATA

OverlapRNA_CD=matrix(NA, nrow=nrow(OverlapRNA), ncol=length(OverlapRNA))
for(i in 1:length(OverlapRNA)) {
	if(as.character(OverlapHeaderRNA["Diet",i])=="CD") {   OverlapRNA_CD[,i]=as.numeric(as.character(OverlapRNA[,i]))	}
} 
rownames(OverlapRNA_CD)=rownames(OverlapRNA)
colnames(OverlapRNA_CD)=colnames(OverlapRNA)

OverlapRNA_HF=matrix(NA, nrow=nrow(OverlapRNA), ncol=length(OverlapRNA))
for(i in 1:length(OverlapRNA)) {
	if(as.character(OverlapHeaderRNA["Diet",i])=="HF") {   OverlapRNA_HF[,i]=as.numeric(as.character(OverlapRNA[,i]))	}
} 
rownames(OverlapRNA_HF)=rownames(OverlapRNA)
colnames(OverlapRNA_HF)=colnames(OverlapRNA)


setwd("/Users/evan.williams/Dropbox/OwnCloud/zy_Shared_Folder/Data/10_EvanAnalysis_Mar2024/")

################ Now we overlap the two datasets and correlate them ##### 
### 
data_in_DNA <- read.table("/Users/evan.williams/Dropbox/OwnCloud/zy_Shared_Folder/Data/10_EvanAnalysis_Mar2024/0_metaDNA_metaphlan4_abundance_all_zy.txt",  header=T, row.names=1, sep="\t")
data_in_RNA <- read.table("/Users/evan.williams/Dropbox/OwnCloud/zy_Shared_Folder/Data/10_EvanAnalysis_Mar2024/0_metaRNA_metaphlan4_abundance_all_zy.txt",  header=T, row.names=1, sep="\t")

# This is just the columns (i.e. samples) that are measured at both nucleotide layers (47 samples)
OverlapDNA=data_in_DNA[,intersect(colnames(data_in_DNA), colnames(data_in_RNA))]
OverlapRNA=data_in_RNA[,intersect(colnames(data_in_RNA), colnames(data_in_DNA))]

OverlapDNA=OverlapDNA[-1,] # Removes the sum total "Bacteria"
OverlapRNA=OverlapRNA[-1,]

rowMeansDNA=rowMeans(OverlapDNA)
rowMeansRNA=rowMeans(OverlapRNA)
isNADNA=NA
isNARNA=NA

for(i in 1:nrow(OverlapDNA)) {
	isNADNA[i]=sum(OverlapDNA[i,]>0.1) # Count how many have reasonable measurements, i.e. at least 0.1%
}

for(i in 1:nrow(OverlapRNA)) {
	isNARNA[i]=sum(OverlapRNA[i,]>0.1) # Count how many have reasonable measurements, i.e. at least 0.1%
}

OverlapDNA[OverlapDNA==0] <- NA # May be better to set 0s to NA
OverlapRNA[OverlapRNA==0] <- NA # May be better to set 0s to NA


## OK, now let's see the correlation between taxa that were reasonably measured at both DNA and RNA levels.
j=1
corrs_rho=NA
corrs_r=NA
pval_rho=NA
pval_r=NA

for(i in 1:nrow(OverlapDNA)) {
	j=which(row.names(OverlapDNA[i,])==row.names(OverlapRNA)) # Find the index overlap
	
	if(rowMeansDNA[i]>0.1 && rowMeansRNA[j]>0.1 && isNADNA[i]>20 && isNARNA[j]>20) {
		corrs_rho[i]=cor(as.numeric(OverlapDNA[i,]),as.numeric(OverlapRNA[j,]), method=c("spearman"), use="complete.obs") 
		corrs_r[i]=cor(as.numeric(OverlapDNA[i,]),as.numeric(OverlapRNA[j,]), method=c("pearson"), use="complete.obs") 
		pval_rho[i]=cor.test(as.numeric(OverlapDNA[i,]),as.numeric(OverlapRNA[j,]),method="spearman")$p.value[1] # Pulls out the p up to e-304
		pval_r[i]=cor.test(as.numeric(OverlapDNA[i,]),as.numeric(OverlapRNA[j,]),method="pearson")$p.value[1] # Pulls out the p up to e-304
	}
	else {
		corrs_rho[i]=NA
		corrs_r[i]=NA
		pval_rho[i]=NA
		pval_r[i]=NA
	}
}

OutTable=cbind(corrs_rho, corrs_r, pval_rho, pval_r)
row.names(OutTable)=row.names(OverlapDNA)

write.table(OutTable,file = "y_CorrTables_MGDNA_MTRNA_AcrossSamples.txt", append = TRUE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names=T, row.names=T, qmethod = c("escape", "double"), fileEncoding = "")

#####################
#### Now let's plot the correlation density at certain taxonomic depths

corrs_rho_phyla=subset(OutTable,!grepl("c__",row.names(OutTable))) # exclude anything with class or below

corrs_rho_classA=subset(OutTable,grepl("c__", row.names(OutTable))) # take everything with class or below
corrs_rho_class=subset(corrs_rho_classA,!grepl("o__",row.names(corrs_rho_classA))) # then exclude anything with order or below

corrs_rho_orderA=subset(OutTable,grepl("o__", row.names(OutTable)))
corrs_rho_order=subset(corrs_rho_orderA,!grepl("f__",row.names(corrs_rho_orderA)))

corrs_rho_familyA=subset(OutTable,grepl("f__", row.names(OutTable)))
corrs_rho_family=subset(corrs_rho_familyA,!grepl("g__",row.names(corrs_rho_familyA)))

corrs_rho_genusA=subset(OutTable,grepl("g__", row.names(OutTable)))
corrs_rho_genus=subset(corrs_rho_genusA,!grepl("s__",row.names(corrs_rho_genusA)))

corrs_rho_speciesA=subset(OutTable,grepl("s__", row.names(OutTable)))
corrs_rho_species=subset(corrs_rho_speciesA,!grepl("t__",row.names(corrs_rho_speciesA)))


plot(density(corrs_rho_phyla[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_phyla[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_phyla[,1]))

plot(density(corrs_rho_class[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_class[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_class[,1]))

plot(density(corrs_rho_order[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_order[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_order[,1]))

plot(density(corrs_rho_family[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_family[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_family[,1]))

plot(density(corrs_rho_genus[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_genus[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_genus[,1]))

plot(density(corrs_rho_species[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_species[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_species[,1]))


################################# FIG S1F SHUFFLED ##########################################
### OK, now we also want to estimate what the false discovery is - so let's shuffle the columns and see how they correlate

OverlapDNA_Rand=sample(OverlapDNA) # Shuffles the columns
OverlapRNA_Rand=sample(OverlapRNA)

## OK, now let's see the correlation between taxa that were reasonably measured at both DNA and RNA levels.
j=1
i=1
corrs_rho_Rand=NA
corrs_r_Rand=NA
pval_rho_Rand=NA
pval_r_Rand=NA

for(i in 1:nrow(OverlapDNA_Rand)) {
	j=which(row.names(OverlapDNA_Rand[i,])==row.names(OverlapRNA_Rand)) # Find the index overlap
	
	if(rowMeansDNA[i]>0.1 && rowMeansRNA[j]>0.1 && isNADNA[i]>20 && isNARNA[j]>20) {
		corrs_rho_Rand[i]=cor(as.numeric(OverlapDNA_Rand[i,]),as.numeric(OverlapRNA_Rand[j,]), method=c("spearman"), use="complete.obs") 
		corrs_r_Rand[i]=cor(as.numeric(OverlapDNA_Rand[i,]),as.numeric(OverlapRNA_Rand[j,]), method=c("pearson"), use="complete.obs") 
		pval_rho_Rand[i]=cor.test(as.numeric(OverlapDNA_Rand[i,]),as.numeric(OverlapRNA_Rand[j,]),method="spearman")$p.value[1] # Pulls out the p up to e-304
		pval_r_Rand[i]=cor.test(as.numeric(OverlapDNA_Rand[i,]),as.numeric(OverlapRNA_Rand[j,]),method="pearson")$p.value[1] # Pulls out the p up to e-304
	}
	else {
		corrs_rho_Rand[i]=NA
		corrs_r_Rand[i]=NA
		pval_rho_Rand[i]=NA
		pval_r_Rand[i]=NA
	}
}

OutTable_Rand=cbind(corrs_rho_Rand, corrs_r_Rand, pval_rho_Rand, pval_r_Rand)
row.names(OutTable_Rand)=row.names(OverlapDNA)

write.table(OutTable_Rand,file = "y_CorrTables_MGDNA_MTRNA_AcrossSamples_RAND.txt", append = TRUE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names=T, row.names=T, qmethod = c("escape", "double"), fileEncoding = "")

#### Now let's plot the FAKE correlation density

plot(density(corrs_rho_Rand, na.rm=T), main="RAND MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 2.2))

### Now let's separate that based on taxonomic level


corrs_rho_rand_phyla=subset(OutTable_Rand,!grepl("c__",row.names(OutTable_Rand))) # exclude anything with class or below

corrs_rho_rand_classA=subset(OutTable_Rand,grepl("c__", row.names(OutTable_Rand))) # take everything with class or below
corrs_rho_rand_class=subset(corrs_rho_rand_classA,!grepl("o__",row.names(corrs_rho_rand_classA))) # then exclude anything with order or below

corrs_rho_rand_orderA=subset(OutTable_Rand,grepl("o__", row.names(OutTable_Rand)))
corrs_rho_rand_order=subset(corrs_rho_rand_orderA,!grepl("f__",row.names(corrs_rho_rand_orderA)))

corrs_rho_rand_familyA=subset(OutTable_Rand,grepl("f__", row.names(OutTable_Rand)))
corrs_rho_rand_family=subset(corrs_rho_rand_familyA,!grepl("g__",row.names(corrs_rho_rand_familyA)))

corrs_rho_rand_genusA=subset(OutTable_Rand,grepl("g__", row.names(OutTable_Rand)))
corrs_rho_rand_genus=subset(corrs_rho_rand_genusA,!grepl("s__",row.names(corrs_rho_rand_genusA)))

corrs_rho_rand_speciesA=subset(OutTable_Rand,grepl("s__", row.names(OutTable_Rand)))
corrs_rho_rand_species=subset(corrs_rho_rand_speciesA,!grepl("t__",row.names(corrs_rho_rand_speciesA)))


plot(density(corrs_rho_rand_phyla[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_rand_phyla[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_rand_phyla[,1]))

plot(density(corrs_rho_rand_class[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_rand_class[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_rand_class[,1]))

plot(density(corrs_rho_rand_order[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_rand_order[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_rand_order[,1]))

plot(density(corrs_rho_rand_family[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_rand_family[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_rand_family[,1]))

plot(density(corrs_rho_rand_genus[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_rand_genus[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_rand_genus[,1]))

plot(density(corrs_rho_rand_species[,1], na.rm=T), main="MG-MT Correlation [Rho]", xlim=c(-1,1), ylim=c(0, 3))
polygon(density(corrs_rho_rand_species[,1], na.rm=T), col="black", border="red")
sum(!is.na(corrs_rho_rand_species[,1]))



## OK, now we have a combined plot showing the random data and the real data. 
plot(density(corrs_rho_phyla[,1], na.rm=T), main="MG-MT Correlation - Phylum Level [Rho]", xlim=c(-1,1), ylim=c(0, 3), las=1,  xaxp=c(-1, 1,4), yaxp=c(0, 3, 3))
polygon(density(corrs_rho_phyla[,1], na.rm=T), col="black", border="red")
lines(density(corrs_rho_rand_phyla[,1], na.rm=T), col="black", lty=2)
sum(!is.na(corrs_rho_phyla[,1]))
sum(!is.na(corrs_rho_rand_phyla[,1]))

#### THIS IS THE EXACT PLOT + DATA EXPORTED
plot(density(corrs_rho_genus[,1], na.rm=T), main="MG-MT Correlation - Genus Level [Rho]", xlim=c(-1,1), ylim=c(0, 3), las=1,  xaxp=c(-1, 1,4), yaxp=c(0, 3, 3))
polygon(density(corrs_rho_genus[,1], na.rm=T), col="blue", border="black")
lines(density(corrs_rho_rand_genus[,1], na.rm=T), col="blue", lty=2)
sum(!is.na(corrs_rho_genus[,1]))
signif(t.test(corrs_rho_genus[,1],corrs_rho_rand_genus[,1], paired=TRUE)$p.value,digits=2) 

GenusCorrs=cbind(corrs_rho_genus[,1], corrs_rho_rand_genus[,1])
row.names(GenusCorrs)=row.names(corrs_rho_genus)
write.table(GenusCorrs,file = "y_CorrTables_MGDNA_MTRNA_AcrossSamples_GenusCorrsOnly.txt", append = TRUE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names=T, row.names=T, qmethod = c("escape", "double"), fileEncoding = "")



























































################################################################################################################
####### Not used in the final figures, but we can show a dot plot of the ratios between the diets
#########################################################################################################

CD=as.numeric(OverlapRNA_CD['k__Bacteria|p__Firmicutes',])/as.numeric(OverlapRNA_CD['k__Bacteria|p__Bacteroidetes',])
HF=as.numeric(OverlapRNA_HF['k__Bacteria|p__Firmicutes',])/as.numeric(OverlapRNA_HF['k__Bacteria|p__Bacteroidetes',])
legendposition="topleft"

groupA=rm.outlier(CD)
groupA=rm.outlier(groupA)
#groupA=groupA[-59]
#groupA=groupA[-64]
#groupA=groupA[-30]
groupB=rm.outlier(HF)
stripchart(groupA, at=0.7, pch=21, cex=1.2, bg="lightblue", method="jitter", xlim=c(0.6, 1.3), las=1, ylim=c(0,25), vertical=T)
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


