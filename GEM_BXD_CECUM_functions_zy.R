BXD.Data.Pre <- function(df, datatype, taxlevel) {
  if (datatype == "mRNA") {
    prodata <- na.omit(df)
    prodata[,2:ncol(prodata)] <- as.numeric(unlist(prodata[,2:ncol(prodata)]))
    return(prodata)
  } else if (datatype == "abundance") {
    prodata <- na.omit(df)
    prodata[,2:ncol(prodata)] <- as.numeric(unlist(prodata[,2:ncol(prodata)]))
    
    tax <- c("phylum", "class", "order", "family", "genus", "species")
    tax_pattern <- c("\\|p__", "\\|c__", "\\|o__", "\\|f__", "\\|g__", "\\|s__", "\\|t__")
    matching_index <- which(tax == taxlevel)
    
    pattern <- tax_pattern[matching_index]
    pattern_exclude <- tax_pattern[matching_index+1]
    
    prodata <- prodata %>%
      subset(grepl(pattern, taxonomy) & !grepl(pattern_exclude, taxonomy))

    pattern_remove <- tax_pattern[matching_index]
    pattern_remove <- gsub('\\\\\\|','',pattern_remove)
    pattern_remove <- paste0(".*\\|(", pattern_remove, ".+)$")

    prodata$taxonomy <- sub(pattern_remove, "\\1", prodata$taxonomy)
    return(prodata)
  }  else if (datatype == "pathabun") {
    prodata <- as.data.frame(df)
    prodata <- na.omit(prodata)
    #prodata <- prodata[-c(1:2),]
    
    prodata <- prodata[!grepl("^UNMAPPED|^UNINTEGRATED", prodata$Pathway), ]
    prodata <- prodata[!grepl("\\|", prodata$Pathway), ]
    prodata[,2:ncol(prodata)] <- as.numeric(unlist(prodata[,2:ncol(prodata)]))
    return(prodata)
  }
}


BXD.ThreeWay.Anova <- function(q, df, meta) {
  All5 <- df
  All <- All5[,-1]
  
  idx <- match(colnames(All), rownames(meta))
  
  Diet <- meta[idx, "Diet"]
  Strain <- meta[idx, "Strain"]
  Age_Categorical <- meta[idx, "AgeC"]
  
  DataRow <- as.numeric(as.character(All[q,]))
  DataRow <- rm.outlier(DataRow, fill=TRUE) 
  
  mainmodel <- aov(DataRow ~ Strain*Diet*Age_Categorical)
  aov_tab <- summary(mainmodel)[[1]]
  ss <- aov_tab[, "Sum Sq"]
  
  GenotypeEffect <- ss[1]	
  DietEffect <- ss[2]		
  AgeEffect <- ss[3]
  GenoDiet <- ss[4]
  StrainAge <- ss[5]
  DietAge <- ss[6]
  GenoDietAge <- ss[7]
  Residual <- ss[8]
  SumSquares <- sum(ss[1:8])
  InteractionEffects <- (GenoDiet + StrainAge + DietAge + GenoDietAge) / SumSquares*100
  
  data.frame(
    ID = as.character(All5[q, 1]),
    GenotypeEffect = GenotypeEffect/SumSquares*100,
    DietEffect = DietEffect/SumSquares*100,
    AgeEffect = AgeEffect/SumSquares*100,
    GenoDiet = GenoDiet/SumSquares*100,
    StrainAge = StrainAge/SumSquares*100,
    DietAge = DietAge/SumSquares*100,
    GenoDietAge = GenoDietAge/SumSquares*100,
    Interaction = InteractionEffects,
    ResidualEffect = Residual/SumSquares*100,
    SumSquares = SumSquares,
    stringsAsFactors = FALSE
  )
}

BXD.Anova.plot <- function(df, meta) {
  res_list <- lapply(seq_len(nrow(df)), function(q) {
    tryCatch(
      BXD.ThreeWay.Anova(q, df, meta),
      error = function(e) { warning(sprintf("Row %d failed: %s", q, e$message)); NULL }
    )
  })
  res_df <- do.call(rbind, res_list)
  res_df <- res_df[complete.cases(res_df), ]
  
  GenotypeEffect <- as.numeric(res_df$GenotypeEffect)
  DietEffect <- as.numeric(res_df$DietEffect)
  AgeEffect <- as.numeric(res_df$AgeEffect)
  Interaction <- as.numeric(res_df$Interaction)
  ResidualEffect <- as.numeric(res_df$ResidualEffect)

  vioplot::vioplot(
    GenotypeEffect, DietEffect, AgeEffect, Interaction, ResidualEffect,
    col = c("deepskyblue","deepskyblue2","deepskyblue4","dodgerblue4","black"),
    ylim = c(0, 100),
    ylab = "Variance Explained [%]",
    xaxt = 'n', yaxt = 'n'
  )
  labels <- c("Strain", "Diet", "Age", "Interact", "Residual")
  axis(1, at = 1:5, labels = labels, lwd = 1, las = 1)
  y_labels <- c("0", "50", "100")
  axis(2, at = c(0, 50, 100), labels = y_labels, lwd = 1, las = 1)
  
  invisible(res_df)
}

BXD.ThreeWay.Fstat = function(q, df, meta) {
  All5 <- df
  All <- All5[,-1]

  idx <- match(colnames(All), rownames(meta))
  
  Diet <- meta[idx, "Diet"]
  Strain <- meta[idx, "Strain"]
  Age_Categorical <- meta[idx, "AgeC"]
  
  DataRow <- as.numeric(as.character(All[q,]))
  DataRow <- rm.outlier(DataRow, fill=TRUE) 
  
  mainmodel <- aov(DataRow ~ Strain*Diet*Age_Categorical)
  aov_tab <- summary(mainmodel)[[1]]
  fv <- aov_tab[, "F value"]
  
  GenotypeEffect <- fv[1]	
  DietEffect <- fv[2]		
  AgeEffect <- fv[3]
  GenoDiet <- fv[4]
  StrainAge <- fv[5]
  DietAge <- fv[6]
  GenoDietAge <- fv[7]
  
  ID <- as.character(All5[q,1])
  
  invisible(data.frame(
    ID = ID,
    GenotypeEffect = fv[1],
    DietEffect = fv[2],
    AgeEffect = fv[3],
    GenoDiet = fv[4],
    StrainAge = fv[5],
    DietAge = fv[6],
    GenoDietAge = fv[7],
    stringsAsFactors = FALSE
  ))
}

BXD.Fstat.plot <- function(df, meta) {
  res_list <- lapply(seq_len(nrow(df)), function(q) {
    tryCatch(BXD.ThreeWay.Fstat(q, df, meta), error = function(e) NULL)
  })
  res_list <- res_list[!sapply(res_list, is.null)]
  if (length(res_list) == 0) stop("All F-stat calculations failed! Please check your input data.")
  
  res_df <- do.call(rbind, res_list)
  res_df <- res_df[complete.cases(res_df), ]
  
  y_max <- max(res_df$GenotypeEffect, res_df$DietEffect, res_df$AgeEffect, na.rm = TRUE)
  
  Fstat_MG_Genotype <- round(res_df$GenotypeEffect, digits = 5)
  Fstat_MG_Diet <- round(res_df$DietEffect, digits = 5)
  Fstat_MG_Age <- round(res_df$AgeEffect, digits = 5)
  
  vioplot(Fstat_MG_Genotype, Fstat_MG_Diet, Fstat_MG_Age,
          col = c("deepskyblue", "deepskyblue2", "deepskyblue4"),
          ylab = "F Statistic",
          ylim = c(0, y_max * 1.05),
          xaxt = 'n')
  labels <- c("Strain", "Diet", "Age")
  axis(1, at = 1:3, labels = labels, lwd = 1, las = 1)
  
  invisible(res_df)
}


BXD.DEG.plot <- function(df, meta, label_n, group = c("Diet", "Age")) {
  group <- match.arg(group)
  rownames(df) <- df[, 1]
  df <- df[, -1]
  idx <- intersect(colnames(df), rownames(meta))
  metadata <- meta[idx, ]
  mRNAdata <- df[, idx]
  
  plot_volcano <- function(res_df, label_n, plot_title) {
    res_df$group <- dplyr::case_when(
      res_df$log2FoldChange > 1 & res_df$padj < 0.05 ~ "Up",
      res_df$log2FoldChange < -1 & res_df$padj < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
    res_df$log10p <- -log10(res_df$pvalue)
    sig_df <- res_df[res_df$group != "NS", ]
    label_genes <- rownames(head(sig_df[order(sig_df$pvalue), ], label_n))
    label_df <- res_df[label_genes, , drop = FALSE]
    p <- ggplot(res_df, aes(x = log2FoldChange, y = log10p, fill = group)) +
      geom_point(shape = 21, color = "black", alpha = ifelse(res_df$group == "NS", 0.2, 0.7), 
                 size = ifelse(res_df$group == "NS", 1.3, 2), stroke = 0.4) +
      scale_fill_manual(values = c("Down" = "#56B4E9", "NS" = "grey80", "Up" = "#E2BE35")) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.4) +
      ggrepel::geom_label_repel(
        data = label_df, aes(label = rownames(label_df)),
        fill = "white", alpha = 0.8, color = "black", fontface = "bold",
        box.padding = 0.35, point.padding = 0.5, label.r = 0.2, label.size = 0.25, size = 3,
        segment.colour = "#4c4b5e"
      ) +
      xlim(-5, 5) +
      ylim(0, min(30, max(res_df$log10p, na.rm = TRUE) * 1.1)) +
      labs(
        title = plot_title,
        x = expression(Log[2]~Fold~Change),
        y = expression(-Log[10]~P~value),
        fill = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_blank()
      )
    print(p)
    res_df
  }
  
  if (group == "Diet") {
    dds <- DESeqDataSetFromMatrix(countData = mRNAdata, colData = metadata, design = ~ Diet + Batch_mRNA)
    dds <- dds[rowSums(counts(dds)) > 1, ]
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("Diet", "HF", "CD"))
    res_df <- as.data.frame(res[complete.cases(res), ])
    res_df_out <- plot_volcano(res_df, label_n, plot_title = "Diet (HF vs CD)")
    return(res_df_out)
  } else if (group == "Age") {
    for (diet_group in c("CD", "HF")) {
      sel <- metadata$Diet == diet_group
      meta_sub <- metadata[sel, ]
      dat_sub <- mRNAdata[, sel]
      dds <- DESeqDataSetFromMatrix(countData = dat_sub, colData = meta_sub, design = ~ AgeC + Batch_mRNA)
      dds <- dds[rowSums(counts(dds)) > 1, ]
      dds <- DESeq(dds)
      res <- results(dds, contrast = c("AgeC", "Old", "Young"))
      res_df <- as.data.frame(res[complete.cases(res), ])
      res_out <- plot_volcano(res_df, label_n, plot_title = paste("Age -", diet_group, "(Old vs Young)"))
      assign(
        x = paste0("degAge", diet_group),
        value = res_out,
        envir = parent.frame()
      )
    }
    invisible(NULL)
  }
}

BXD.ORA.plot <- function(df, top_n) {
  df$SYMBOL <- rownames(df) 
  
  gene_df <- bitr(df$SYMBOL,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
  
  res_annot <- merge(df, gene_df, by = "SYMBOL")
  res_annot <- res_annot[order(res_annot$ENTREZID), ]
  
  deg <- res_annot[res_annot$log2FoldChange > 1 & res_annot$padj < 0.05 |
                     res_annot$log2FoldChange < -1 & res_annot$padj < 0.05, ]
  deg_genes <- unique(deg$ENTREZID)
  
  set.seed(12345)
  kegg_ora <- enrichKEGG(gene = deg_genes,
                         organism = "mmu",
                         pvalueCutoff = 0.1,
                         qvalueCutoff = 0.1)
  kegg_ora <- setReadable(kegg_ora, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  kegg_ora@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_ora@result$Description)
  
  top_terms <- kegg_ora@result[order(kegg_ora@result$qvalue), ][1:top_n, ]
  top_terms$Description <- factor(top_terms$Description,
                                  levels = rev(top_terms$Description))  
  
  ggplot(top_terms, aes(x = RichFactor,
                        y = Description,
                        size = Count,
                        fill = qvalue)) +
    geom_point(shape = 21,
               color = "black",
               aes(fill = qvalue),
               stroke = 0.3) +
    coord_cartesian(xlim = c(min(top_terms$RichFactor) - 0.02,
                             max(top_terms$RichFactor) + 0.02)) +
    theme_minimal(base_size = 14) +
    labs(x = "Rich Factor", y = "", size = "Gene Number") +
    theme(axis.text.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.title = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          panel.grid = element_line(color = "grey90"),
          axis.ticks = element_line(color = "black", size = 0.3),
          axis.ticks.length = unit(0.1, "cm"))
}


BXD.Ratio.plot <- function(df, meta, tax1, tax2) {
  df_T1 <- df[grep(paste0("^", tax1), df$taxonomy), ]
  df_T2 <- df[grep(paste0("^", tax2), df$taxonomy), ]
  
  T1T2_Ratio <- df_T1[, 2:ncol(df_T1)] / df_T2[, 2:ncol(df_T2)]
  t1t2name <- paste0(tax1, "/", tax2)
  T1T2_Ratio <- cbind(taxonomy = t1t2name, T1T2_Ratio)
  T1T2R <- rbind(df, T1T2_Ratio)
  
  meta$sample <- rownames(meta)
  
  T1T2R_long <- 
    pivot_longer(T1T2R, 
                 cols = -contains("taxonomy"), # columns argument, required
                 names_to = "sample",
                 values_to = "abundance")
  
  T1T2R_long <- left_join(T1T2R_long,
                          meta %>%
                            select(sample, Diet, TIMEonDIET, Age_group, Age),
                          by = "sample")
  T1T2R_long$agediet <- paste(T1T2R_long$Age_group, T1T2R_long$Diet, sep = "_")
  
  T1T2R <- T1T2R_long %>% subset(taxonomy == t1t2name)
  
  p1 <- ggplot(T1T2R, aes(x = Diet, y = abundance, fill = Diet)) +
    geom_jitter(aes(col = Diet),
                alpha = 0.7, 
                size = 4,
                position = position_jitterdodge(jitter.width = 0.35,
                                                jitter.height = 0,
                                                dodge.width = 0.8)) +
    geom_boxplot(width=0.5,
                 alpha = 0.4,
                 position = position_dodge(0.1),
                 size = 0.4, outlier.colour = NA, fill = 'white') +
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_y_continuous(limits = c(0,25), breaks = c(0,5,10,20)) +
    scale_color_manual(values = c('#56B4E9','#E2BE35')) +
    stat_compare_means(aes(group= Diet),
                       comparisons = list(c("CD", "HF")),
                       method = "t.test",
                       label="p.value",   
                       label.x = 1.4,
                       label.y = 0.95)+
    labs(x = "", y = "Ratio") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16), 
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)) 
  
  p2 <- ggplot(T1T2R, aes(x = agediet, y = abundance, fill = Diet)) +
    geom_jitter(aes(col = Diet),
                alpha = 0.7, 
                size = 4,
                position = position_jitterdodge(jitter.width = 0.35,
                                                jitter.height = 0,
                                                dodge.width = 0.8)) +
    geom_boxplot(width=0.5,
                 alpha = 0.4,
                 position = position_dodge(0.1),
                 size = 0.4, outlier.colour = NA, fill = 'white') +
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_y_continuous(limits = c(0,25), breaks = c(0,5,10,20)) +
    scale_color_manual(values = c('#56B4E9','#E2BE35')) +
    stat_compare_means(aes(group= agediet),
                       comparisons = list(c("AG1_CD", "AG1_HF"), c("AG2_CD", "AG2_HF"), c("AG3_CD", "AG3_HF"), c("AG4_CD", "AG4_HF")),
                       method = "t.test",
                       label="p.value",   
                       label.x = 1.4,
                       label.y = 0.95)+
    labs(x = "", y = "Ratio") +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 16), 
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)) 
  
  print(p1)
  print(p2)
}

BXD.Stacked.data <- function(df_sum, meta, classification, top_taxa, dfType) {
  topN_common <- 10
  topN_unique <- 15
  group <- "Diet"
  
  df_other <- colSums(df_sum[(!df_sum[[classification]] %in% top_taxa), 2:ncol(df_sum)])
  df_other <- as.data.frame(df_other)
  df_other <- t(df_other)
  df_other <- cbind("Others", df_other)
  colnames(df_other)[1] <- classification
  df_sort <- df_sum[(df_sum[[classification]] %in% top_taxa), ]
  df_sort <- rbind(df_sort, df_other)
  
  sub_design <- meta[order(meta$Diet,meta$Age,meta$Strain),]
  df_idx <- rownames(sub_design) %in% colnames(df_sort) 
  df_sub_design <- sub_design[df_idx,]
  df_list <- rownames(df_sub_design)
  df_sample_order <- c(df_list)
  df_sub_design$Sample <- rownames(df_sub_design)
  
  df_sub_design[,group] <- as.factor(df_sub_design[,group])
  df_g1 <- levels(df_sub_design[,group])[1]
  df_g2 <- levels(df_sub_design[,group])[2]
  df_g1_list <- rownames(df_sub_design[df_sub_design[,group] %in% df_g1,])
  df_g2_list <- rownames(df_sub_design[df_sub_design[,group] %in% df_g2,])
  df_sort_g1 <- df_sort[,df_g1_list]
  df_sort_g2 <- df_sort[,df_g2_list]
  rownames(df_sort_g1) <- df_sort[,1]
  rownames(df_sort_g2) <- df_sort[,1]
  df_sort_g1[1:ncol(df_sort_g1)] <- as.numeric(unlist(df_sort_g1[1:ncol(df_sort_g1)]))
  df_sort_g2[1:ncol(df_sort_g2)] <- as.numeric(unlist(df_sort_g2[1:ncol(df_sort_g2)]))
  df_p_val <- NULL
  for (taxon in rownames(df_sort_g1)) {
    p <- t.test(df_sort_g1[taxon,], df_sort_g2[taxon,])$p.value
    df_p_val <- rbind(df_p_val, c(taxon, p))
  }
  colnames(df_p_val) <- c("classification", "p")
  df_p_val <- as.data.frame(df_p_val)
  df_p_val$p <- as.numeric(format(df_p_val$p, scientific = FALSE))
  df_p_val <- df_p_val[order(df_p_val$p),]
  others_row <- which(df_p_val$classification == "Others")
  df_p_val <- rbind(df_p_val, df_p_val[others_row,])
  df_p_val <- df_p_val[-others_row,]
  df_tax_list <- rev(df_p_val[,1])
  df_tax_order <- c(df_tax_list)
  
  df_long <- 
    pivot_longer(df_sort, 
                 cols = !all_of(classification),
                 names_to = "Sample",
                 values_to = "RA")
  names(df_long) <- c("classification", "Sample", "RA")
  df_long[,3] <- as.numeric(unlist(df_long[,3]))
  df_long <- left_join(df_long, df_sub_design %>% select(Sample, Diet), by = "Sample")
  df_long$Diet <- paste0(dfType, df_long$Diet)
  
  df_long <- df_long %>% 
    group_by(Sample) %>% 
    mutate(total = sum(RA)) %>%
    ungroup() %>%
    mutate(relative_abundance = RA / total)
  df_long$classification <- factor(df_long$classification,
                              levels = df_tax_order)
  df_long$Sample <- factor(df_long$Sample,
                            levels = df_sample_order)

  return(df_long)
}

BXD.Stacked.plot <- function(df1, df2, meta, classification, dfType1, dfType2) {
  topN_common <- 10
  topN_unique <- 15
  group <- "Diet"
  
  df1_sum <- df1[(order(-rowSums(df1[,2:ncol(df1)]))),] 
  df2_sum <- df2[(order(-rowSums(df2[,2:ncol(df2)]))),] 
  top_common_taxa <- head(intersect(df1[[classification]], df2[[classification]]), topN_common)
  df1_top_taxa <- setdiff(head(df1_sum[[classification]], topN_unique), top_common_taxa)
  df2_top_taxa <- setdiff(head(df2_sum[[classification]], topN_unique), top_common_taxa)
  top_taxa <- c(top_common_taxa, df1_top_taxa, df2_top_taxa)
  topN <- length(top_taxa) + 1

  df1_long <- BXD.Stacked.data(df1_sum, meta, classification, top_taxa, dfType1)
  df2_long <- BXD.Stacked.data(df2_sum, meta, classification, top_taxa, dfType2)
  
  sum_sort_long <- rbind(df1_long, df2_long)
  
  ggplot(sum_sort_long, aes(x = Sample,
                            y = relative_abundance,
                            fill = classification,
                            alluvium = classification
  )) +
    geom_bar(position = "fill",
             stat = "identity",
             color = "black",
             size = 0.05,
             width = 0.7) +
    scale_y_continuous(expand = c(0, 0), 
                       labels = scales::percent_format()) +
    labs(x = "", y = "Relative Abundance [%]") +
    facet_grid(~Diet, scales = "free_x", space = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, ncol = 3))
}

