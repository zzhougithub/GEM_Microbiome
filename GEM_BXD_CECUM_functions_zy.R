BXD.Data.Pre <- function(df, datatype, taxlevel) {
  if (datatype == "mRNA") {
    prodata <- na.omit(df)
    prodata[,2:ncol(prodata)] <- as.numeric(unlist(prodata[,2:ncol(prodata)]))
    return(prodata)
  } else if (datatype %in% c("MG", "MT")) {
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


BXD.DEG.plot <- function(df, meta, label_n = 20, group = c("Diet", "Age")) {
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

