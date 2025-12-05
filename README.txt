Getting started

This is the code and supplemental data files for the revised submission to Microbiome. The original 2025 paper on bioRxiv is: https://www.biorxiv.org/content/10.1101/2025.07.31.667922v1.abstract

This README will be updated with the final published details of the manuscript, while all other files will remain as-is from the final revised submission.

** Manuscript Name

Original: Integrative Analysis of the Mouse Cecal Microbiome Across Diet, Age, and Metabolic State in the Diverse BXD Population
Revised: Integrative Analysis of the Mouse Cecal Microbiome Across Diet, Age, and Weight in the Diverse BXD Population

——————————————————————————————
** Manuscript Description

The gut microbiota both adapts to and shapes the host’s metabolic state through affecting circulating metabolites and consequent gene regulatory networks, resulting in systemic influences in diverse organs via connections such as the gut–liver axis. Numerous variables such as diet, age, and host genetics modulate the composition of the gut microbiome, but their interactions and specific associative and mechanistic links to host molecular phenotypes remain incompletely unannotated and, consequently, poorly understood. Integrated multi-omics approaches in genetically diverse populations offer an opportunity to dissect these interactions and identify predictive microbial signatures for host phenotypes, whether systemic such as body weight, or specific such as the molecular activity of metabolic pathways in gut tissue.

We profiled the cecal metagenome, metatranscriptome, and host transcriptome from 232 mice belonging to 175 distinct cohorts according to a low fat chow diet (CD) or a high-fat diet (HF), four adult ages (between roughly 180 to 730 days of age), and 43 distinct genotypes (all different, inbred BXD strains). Genetics and diet exerted the strongest influence on microbiota abundance and activity, followed distantly by age. HF feeding significantly reduced diversity across all ages and all genotypes, altering >300 species. Machine learning models based on microbial profiles reliably predicted body weight within dietary group (AUC = 0.87 for CD, 0.83 for HF) and chronological age (AUC = 0.84), with model performance rising to 0.95 when integrating top microbial features with liver proteomics. Network analyses of expression data revealed links between genes, pathways, and specific microbes, including a negative association between cecal Ido1 expression and short-chain fatty acid (SCFA)-producing Lachnospiraceae, suggesting dietary fat may modulate host tryptophan metabolism through microbiota shifts.

Broadly, it is necessary to advance the community’s understanding of the microbiome-host interactions to identify compact panels of biomarkers that predict outcomes such as body weight. The large contribution of genetics, diet, and other environmental parameters to variation in the microbiota abundance and activity cause challenges for creating generalizable panels of highly-predictive biomarkers. In our cross-sectional multi-omics data from aging mouse, we identify highly predictive and robust microbial biomarkers of aging and obesity. Yet, this is under ideal circumstances, as the models were trained on mice with similar genotypes, diets, and ages. Challenges will certainly arise when applying this approach to outbred or clinical populations, yet we here establish that microbiome data are sufficiently reliable to provide high predictive capacities even when major factors like genetic background, diet, and age are variable.


————————————————————————————————————————————————————————————

** File Description

- data - This folder contains the supplemental tables for the manuscript, as follows:
	• S1_metadata_cecalsamples.xlsx - A two sheet XLSX file. The first sheet is all ceca that were collected and frozen in the aging colony. The second sheet is all ceca selected for MG, MT, and/or cecal mRNA analysis in this paper.
	• S2_mRNA_cecum_raw_counts_zy.txt - This is the raw counts from the cecal mRNA samples, from the normal mouse tissue.
	• S3_mRNA_metaphlan4_abundance_zy.txt - The same file as above but adjusted by TPM and scaled to log2.
	• S4_metaDNA_metaphlan4_abundance_zy.txt - This is the results at all taxonomic ranks of relative abundances in MG data across all samples as from the metaphlan4 pipeline. Note that to perform analyses, you should subset and only take a particular taxonomic rank (e.g. all phyla, or all classes - do not use this entire table at once). Also note that technical replicate samples (those with suffix _r or _n) should be removed prior to running analyses on machine learning. 
	• S5_metaRNA_metaphlan4_abundance_zy.txt - This is the results at all taxonomic ranks of relative abundances in MT data across all samples as from the metaphlan4 pipeline. Note that to perform analyses, you should subset and only take a particular taxonomic rank (e.g. all phyla, or all classes - do not use this entire table at once).
	• S6_metaDNA_humann3.8_pathabundance_relab-zy.txt - The functional pathway 'abundances' calculated with humann3 using MG data.
	• S7_metaRNA_humann3.8_pathabundance_relab-zy.txt - The functional pathway 'activities' calculated with humann3 using MT data.
	• S8_metaRNA_SGB2GTDB_zy.txt - Lookup table between SGB numerical identifiers and the alternative GTDB numerical identifiers. We expect that many of these poorly-annotated species will have improved annotation in coming years as microbiome data becomes more complete.
	• S9_metaDNA_kraken2_abundance_zy.txt - This is equivalent to Table S4, but the output from the kraken2 pipeline instead of metaphlan. Taxonomic identifications still can vary quite significantly between pipelines, especially at the lowest taxonomic rank (i.e. species).
	
	
- 1_rawdata_MG&MT - This is the basic shell code to run the kneaddata, metaphlan, and humann pipelines from the MG and MT fastq files to generate the supplemental tables.
- 1_rawdata_Transcriptome.R - This is the basic R code to run Rsubread to perform mRNA count analysis from the cecal transcriptome fastq files to generate the supplemental tables. 
- 2_GEM_BXD_CECUM_2025_zy.R - This file contains the main R code used to directly make R-derived figures in the paper. Note that some figures are hand-drawn (e.g. Fig 1A, Fig 1B) and other figures are generated from GUI-based software (e.g. qPCR analysis in Figs S1B-C-D, or GSEA figures from Fig S2 C and E).
- 2_GEM_BXD_CECUM_functions_zy.R - This file contains some support functions that are not part of standard R libraries which are called in the above file, but which do not need to be edited or modified.
- 2_MGMTCorrelations_FigS1F.R - This is the R code that was used to check the Spearman and Pearson correlations between the MG and MT abundances for the 48 overlapping samples, for the 38 most abundant microbes. Plotted in FigS1F, showing that generally MG and MT data agree in terms of abundance calculations of microbes in the gut microbiome - but not always.
- 2_VariancesExplained_Code.r - Code used to calculate variance explained according to key study variables (e.g. technical replicates, extraction replicates, biological replicates, all samples...)
- 3_GEM_BXD_CECUM_ML_AGE_zy.py - Python scripts for developing and testing models and performing feature selection and calculating feature importance in Figs 4 and 5 and their associated supplemental figures.
- 3_GEM_BXD_CECUM_ML_utils_zy.py - Support scripts loaded as part of the above python script; no modifications should be necessary here.
- 4_Regression_Correlation_Lachnospiraceae_M18_1_Ido1.Rmd - Regression analysis used to examine the relationship between Ido1 and a particular species of Lachnospiraceae across and within diet, and with or without mediation. While diet has a strong effect on both variables, there is a link between Ido1 and that Lachnospiraceae independently of diet, particularly in CD cohorts, and is observed at both the MG and MT level.
- LICENSE.txt - This is the license file that allows full use of these scripts and data.
- README.md - This readme file.

Also note that analyses for Figure 5 rely on using liver transcriptome and liver proteome generated in the same animals but in a prior publication of ours. These data are accessible in Dataset S1 in the manuscript "Multiomic profiling of the liver across diets and age in a diverse mouse population" (Williams et al., Cell Systems, 2022). This paper is open access. A direct link to that file is here: https://ars.els-cdn.com/content/image/1-s2.0-S2405471221003446-mmc8.zip (44 MB)

Most code by Ziyun Zhou, some contributions by Evan Williams and Arianna Lamanna. 
Other direct work on data in the manuscript, including data tables, Arianna Lamanna, Rashi Halder, Emeline Pansart, Besma Boussoufa, Thamila Kerkour, and Evan G. Williams. Other authors (Paul Wilmes, Shaman Narayanasamy) did not work on code or data generation, but worked on text and interpretation. We would like to give appreciation to all of those who have contributed to the project, such as the initial generation of the aging cohort and the tissue collection, particularly Robert Williams and Suheeta Roy at UTHSC Memphis.

Project is complete, and this code will be frozen to pair with the revised submission of this manuscript in December 2025. 
