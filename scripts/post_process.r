#title: Post-process
#author: Delaney Smith
#date: '2025-6-25'
#function: fit combined (DNA + DNAm) model per gene, collect R2, feature counts and weights for all models, test fold-wise R2 differences
#output: combined model files in combined_model_output, summaries in Output 4 directory (pgx_ prefix for PGx genes, none for the rest), all_results_gene_stats.csv and pgx_all_results_gene_stats.csv in working directory

#load libraries
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(tidyverse)
library(minfi)
library(data.table)
library(glmnet)
library(stringr)

##################### set up ################################

#file paths
base_dir <- "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4"
dna_dir <- "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/DNA_model_output"
dnam_dir <- "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/DNAm_model_output"
out_dir <- "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/combined_model_output"

pgx_genes <- c("CYP2B6", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A5", "DPYD", "F5", "HLA-B", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1")

#genes with DNA, DNAm and combined output (defines the non-PGx gene set)
dna_files  <- list.files(dna_dir,  pattern = "_DNA\\.RData$",   full.names = FALSE)
dnam_files <- list.files(dnam_dir, pattern = "_DNAm\\.RData$",  full.names = FALSE)
comb_files <- list.files(out_dir,  pattern = "_comb\\.RData$",  full.names = FALSE)
dna_genes  <- sub("_DNA\\.RData$", "", dna_files)
dnam_genes <- sub("_DNAm\\.RData$", "", dnam_files)
comb_genes <- sub("_comb\\.RData$", "", comb_files)
common_genes <- Reduce(intersect, list(dna_genes, dnam_genes, comb_genes))

#PGx genes, or all other modeled genes (plus HLA-B, which is in both sets)
in_gene_set <- function(gene_name, gene_set) {
  if (gene_set == "pgx") {
    gene_name %in% pgx_genes
  } else {
    gene_name %in% common_genes && (!gene_name %in% pgx_genes || gene_name == "HLA-B")
  }
}

#run the full pipeline for the non-PGx genes (no file prefix) and the PGx genes (pgx_ prefix)
for (gene_set in c("normal", "pgx")) {

  prefix <- if (gene_set == "pgx") "pgx_" else ""
  cat("Gene set:", gene_set, "\n")

  #load expression (y) data
  if (gene_set == "pgx") {
    load("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/rpkm_pgx.RData")
    rpkm_se <- rpkm_pgx
  } else {
    load("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/rpkm_pgkb.RData")
  }

  ##################### make combined model ################################

  #per gene: stack out-of-fold DNA and DNAm predictions in an elastic net
  #saves <gene>_comb.RData (R2, cvfit, null flags) to out_dir
  dna_files  <- list.files(dna_dir,  pattern = "_DNA\\.RData$",   full.names = TRUE)

  for (dna_file in dna_files) {

    gene_name <- sub("_DNA\\.RData$", "", basename(dna_file))

    if (!in_gene_set(gene_name, gene_set)) next

    dnam_file <- file.path(dnam_dir, paste0(gene_name, "_DNAm.RData"))

    load(dna_file)
    d_cvfit <- DNA_results$cvfit
    load(dnam_file)
    m_cvfit <- DNAm_results$cvfit

    #lambda.min index and whether the model has any non-zero coefficients
    id_dna  <- which.min(abs(d_cvfit$glmnet.fit$lambda - d_cvfit$lambda.min))
    id_dnam <- which.min(abs(m_cvfit$glmnet.fit$lambda - m_cvfit$lambda.min))

    dna_null  <- sum(coef(d_cvfit$glmnet.fit, s = d_cvfit$lambda[id_dna])[-1] != 0) == 0
    dnam_null <- sum(coef(m_cvfit$glmnet.fit, s = m_cvfit$lambda[id_dnam])[-1] != 0) == 0

    if (dna_null || dnam_null)
      message(gene_name, ": DNA null = ", dna_null, ", DNAm null = ", dnam_null)

    #skip combined model if either model is null (R2 = NA)
    if (dna_null || dnam_null) {
      combined_results <- list(
        EnsemblId = DNAm_results$EnsemblId,
        R2 = NA,
        cvfit = NULL,
        DNA_null = dna_null,
        DNAm_null = dnam_null
      )
      save(combined_results, file = file.path(out_dir, paste0(gene_name, "_comb.RData")))
      next
    }

    #out-of-fold predictions aligned on shared sample IDs
    pred_dna  <- as.numeric(d_cvfit$fit.preval[, id_dna])
    pred_dnam <- as.numeric(m_cvfit$fit.preval[, id_dnam])

    ids_dna  <- rownames(d_cvfit$fit.preval)
    ids_dnam <- rownames(m_cvfit$fit.preval)
    common_ids <- intersect(ids_dna, ids_dnam)

    X <- cbind(
      DNA  = pred_dna[match(common_ids, ids_dna)],
      DNAm = pred_dnam[match(common_ids, ids_dnam)]
    )

    #expression of this gene for the shared samples
    y_gene <- rpkm_se[rownames(rpkm_se) == DNAm_results$EnsemblId, common_ids, drop = FALSE]
    y <- as.numeric(t(as.matrix(assay(y_gene))))

    #reuse the DNA model folds
    foldid <- d_cvfit$foldid[match(common_ids, ids_dna)]

    #combined model and its out-of-fold R2
    meta_cv <- cv.glmnet(as.matrix(X), y, alpha = 0.5, keep = TRUE, foldid = foldid)

    id_meta <- which.min(abs(meta_cv$lambda - meta_cv$lambda.min))
    preds_meta <- as.numeric(meta_cv$fit.preval[, id_meta])
    R2_meta <- cor(preds_meta, y)^2

    combined_results <- list(
      EnsemblId = DNAm_results$EnsemblId,
      R2 = R2_meta,
      cvfit = meta_cv,
      DNA_null = dna_null,
      DNAm_null = dnam_null
    )

    save(combined_results, file = file.path(out_dir, paste0(gene_name, "_comb.RData")))
    gc()
  }

  ############# process combined results #############

  #per model type: collect gene, EnsemblId and R2 into a data frame
  #saves <prefix>summary_results_<model>.RData to base_dir
  model_types <- c("DNAm", "DNA", "combined")

  for (model in model_types) {

    cat("Processing model type:", model, "\n")

    #set paths
    data_dir <- file.path(base_dir, paste0(model, "_model_output"))
    files <- list.files(data_dir, pattern = "\\.RData$", full.names = TRUE)
    file_names <- basename(files)
    gene_names <- sub("_.*", "", file_names)

    #initialize list
    results_list <- list()

    #load and extract results
    for (i in seq_along(files)) {

      f <- files[i]
      gene_name <- gene_names[i]
      if (!in_gene_set(gene_name, gene_set)) next

      load(f)
      varname <- paste0(model, "_results")

      if (exists(varname)) {
        obj <- get(varname)
        results_list[[length(results_list) + 1]] <- list(
          Gene = gene_name,
          EnsemblId = obj$EnsemblId,
          R2 = obj$R2
        )
        rm(list = varname)
      } else {
        warning(paste("Results object", varname, "not found in", f))
      }

      gc()
    }

    #combine into data frame
    df <- data.frame(
      Gene = sapply(results_list, `[[`, "Gene"),
      EnsemblId = sapply(results_list, `[[`, "EnsemblId"),
      R2 = sapply(results_list, `[[`, "R2"),
      stringsAsFactors = FALSE
    )

    #check for and drop nulls
    na_idx <- is.na(df$R2)
    if (any(na_idx)) {
      cat("Removing genes with NA R2:\n")
      cat(paste(df$Gene[na_idx], collapse = ", "), "\n")
    }

    df <- df[!na_idx, , drop = FALSE]

    #save
    save_path <- file.path(base_dir, paste0(prefix, "summary_results_", model, ".RData"))
    save(df, file = save_path)
    cat("Saved summary to:", save_path, "\n")

  }

  ############## combine all modalities and reformat #############

  #start here if previous steps already done to save time
  load(file = file.path(base_dir, paste0(prefix, "summary_results_DNAm.RData")))
  DNAm_df <- df
  load(file = file.path(base_dir, paste0(prefix, "summary_results_DNA.RData")))
  DNA_df <- df
  load(file = file.path(base_dir, paste0(prefix, "summary_results_combined.RData")))
  comb_df <- df

  #reformat
  DNAm_df$group <- "DNAm"
  DNA_df$group <- "Genotype"
  comb_df$group <- "Combined"
  results_df <- bind_rows(DNAm_df, DNA_df, comb_df)

  #rename R2 columns
  colnames(DNAm_df)[colnames(DNAm_df) == "R2"] <- "R2_DNAm"
  colnames(DNA_df)[colnames(DNA_df) == "R2"] <- "R2_DNA"
  colnames(comb_df)[colnames(comb_df) == "R2"] <- "R2_comb"

  #merge into one wide data frame, one row per gene
  merged_df <- DNAm_df %>%
    inner_join(DNA_df,  by = c("Gene", "EnsemblId")) %>%
    inner_join(comb_df, by = c("Gene", "EnsemblId")) %>%
    select(Gene, EnsemblId, R2_DNAm, R2_DNA, R2_comb)

  save(merged_df, file = file.path(base_dir, paste0(prefix, "all_results_1.RData")))

  ############# add number of features for each model #############

  #per model (DNAm, DNA): non-zero coefficient count and mean absolute effect at lambda.min
  #produces DNA_counts, DNAm_counts and summary_list (saved as <prefix>model_feat_sum_stats.RData)
  model_types <- c("DNAm", "DNA")
  DNA_counts  <- data.frame(Gene = character(), n_features = numeric(), stringsAsFactors = FALSE)
  DNAm_counts <- data.frame(Gene = character(), n_features = numeric(), stringsAsFactors = FALSE)

  summary_list <- list()

  for (model in model_types) {
    cat("\nProcessing model type:", model, "\n")

    dir_path <- file.path(base_dir, paste0(model, "_model_output"))
    pattern  <- paste0(model, "\\.RData$")
    rdata_files <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)

    counts      <- c()
    r2_values   <- c()
    abs_effects <- c()
    genes       <- c()

    for (i in seq_along(rdata_files)) {
      gene_name <- sub(paste0("_", model, "\\.RData$"), "", basename(rdata_files[i]))

      if (!in_gene_set(gene_name, gene_set)) next

      load(rdata_files[i])

      varname <- paste0(model, "_results")
      if (exists(varname)) {
        obj <- get(varname)
      } else {
        warning(paste("Results object", varname, "not found in", rdata_files[i]))
      }

      cvfit <- obj$cvfit
      lambda_idx  <- which.min(abs(cvfit$glmnet.fit$lambda - cvfit$lambda.min))
      coef_sparse <- cvfit$glmnet.fit$beta[, lambda_idx]
      nonzero_coef <- coef_sparse[coef_sparse != 0]

      #deal with null models (shouldn't be happening now)
      if (length(nonzero_coef) == 0) {
        cat("Dropped null model (all-zero coefficients):", gene_name, "\n")
        counts      <- c(counts, 0)
        r2_values   <- c(r2_values, obj$R2)
        abs_effects <- c(abs_effects, NA_real_)
        genes       <- c(genes, gene_name)
        rm(list = varname)
        gc()
        next
      }

      coef_values <- as.numeric(nonzero_coef)

      #append results
      counts      <- c(counts, length(nonzero_coef))
      r2_values   <- c(r2_values, obj$R2)
      abs_effects <- c(abs_effects, mean(abs(coef_values)))
      genes       <- c(genes, gene_name)

      rm(list = varname)
      gc()
    }

    #per-gene feature counts for this model
    if (length(genes) > 0) {
      model_counts <- data.frame(
        Gene = genes,
        n_features = counts,
        stringsAsFactors = FALSE
      )

      if (model == "DNA") {
        DNA_counts <- model_counts
      } else if (model == "DNAm") {
        DNAm_counts <- model_counts
      }
    }

    #summary statistics
    stats <- list(
      model = model,
      mean_nonzero       = mean(counts, na.rm = TRUE),
      range_nonzero      = range(counts, na.rm = TRUE),
      cor_nonzero_r2     = cor(counts, r2_values, use = "complete.obs"),
      mean_abs_effect    = mean(abs_effects, na.rm = TRUE),
      cor_abs_effect_r2  = cor(abs_effects, r2_values, use = "complete.obs")
    )

    summary_list[[model]] <- stats

    cat(sprintf(
      "%s model: mean # of features = %.1f, range = [%d, %d], cor(feats,R2)=%.3f\n",
      model,
      stats$mean_nonzero,
      stats$range_nonzero[1],
      stats$range_nonzero[2],
      stats$cor_nonzero_r2
    ))
    cat(sprintf(
      "mean abs(effect)=%.5f, cor(abs_effect,R2)=%.3f\n",
      stats$mean_abs_effect,
      stats$cor_abs_effect_r2
    ))
  }

  save(summary_list, file = file.path(base_dir, paste0(prefix, "model_feat_sum_stats.RData")))

  #add feature counts to merged_df
  DNA_counts_tmp  <- dplyr::rename(DNA_counts,  n_DNA_features  = n_features)
  DNAm_counts_tmp <- dplyr::rename(DNAm_counts, n_DNAm_features = n_features)

  merged_df <- merged_df %>%
    left_join(DNA_counts_tmp,  by = "Gene") %>%
    left_join(DNAm_counts_tmp, by = "Gene") %>%
    mutate(
      n_DNA_features  = ifelse(is.na(n_DNA_features),  0, n_DNA_features),
      n_DNAm_features = ifelse(is.na(n_DNAm_features), 0, n_DNAm_features)
    )

  save(merged_df, file = file.path(base_dir, paste0(prefix, "all_results_2.RData")))

  ############# add combined model weights #############

  #extract DNA and DNAm weights from each combined model, then relative weights
  #adds w_dna, w_dnam, prop_dna and prop_dnam to merged_df
  files <- list.files("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/combined_model_output/", pattern = "comb\\.RData$", full.names = TRUE)

  meta_weights_list <- lapply(files, function(f) {
    gene_name <- sub("_comb\\.RData$", "", basename(f))

    if (!in_gene_set(gene_name, gene_set)) return(NULL)

    combined_results <- get(load(f))
    cvfit <- combined_results$cvfit

    #null combined models have no cvfit and no weights
    if (is.null(cvfit)) return(NULL)

    #lambda.min index
    id <- which(cvfit$lambda == cvfit$lambda.min)

    #extract coefficients (intercept + two modalities)
    coefs <- as.numeric(coef(cvfit, s = cvfit$lambda.min))
    names(coefs) <- rownames(coef(cvfit))

    #coefs = (Intercept, pred_dna, pred_dnam)
    data.frame(
      EnsemblId = combined_results$EnsemblId,
      w_dna = coefs["DNA"],
      w_dnam = coefs["DNAm"]
    )
  })

  meta_weights_list <- Filter(Negate(is.null), meta_weights_list)
  meta_weights_df <- do.call(rbind, meta_weights_list)

  merged_df <- merged_df %>%
    left_join(
      meta_weights_df %>%
        select(EnsemblId, w_dna, w_dnam),
      by = "EnsemblId"
    )

  #relative weight of each modality
  merged_df <- merged_df %>%
    mutate(
      w_total = abs(w_dna) + abs(w_dnam),
      prop_dna = abs(w_dna) / w_total,
      prop_dnam = abs(w_dnam) / w_total
    ) %>%
    select(-w_total)

  save(merged_df, file = file.path(base_dir, paste0(prefix, "all_results_final.RData")))

  ############# gene-wise significance testing for performance #############

  #fold-wise R2 of one model for one gene
  #group: "Genotype", "DNAm" or "meta"; returns data frame of fold, R2, group, Symbol
  extract_fold_r2 <- function(rpkm_se, group, symbol) {
    if (group == "Genotype") {
      file_path <- sprintf("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/DNA_model_output/%s_DNA.RData", symbol)
      results <- get(load(file_path))
      cvfit <- results$cvfit
      ids_dna  <- rownames(cvfit$fit.preval)
    } else if (group == "DNAm") {
      file_path <- sprintf("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/DNAm_model_output/%s_DNAm.RData", symbol)
      results <- get(load(file_path))
      cvfit <- results$cvfit
      ids_dna  <- rownames(cvfit$fit.preval)
    } else if (group == "meta") {
      file_path <- sprintf("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/combined_model_output/%s_comb.RData", symbol)
      results <- get(load(file_path))
      cvfit <- results$cvfit

      name_path <- sprintf("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/DNA_model_output/%s_DNA.RData", symbol)
      names <- get(load(name_path))
      n <- names$cvfit
      ids_dna  <- rownames(n$fit.preval)
    } else {
      stop("Unsupported group: must be 'Genotype', 'DNAm', or 'meta'")
    }

    common_ids <- intersect(ids_dna, colnames(rpkm_se))
    ids_dna <- common_ids

    id <- which(cvfit$lambda == cvfit$lambda.min)
    preds <- cvfit$fit.preval[, id]
    foldid <- cvfit$foldid
    y_gene <- rpkm_se[rownames(rpkm_se) == results$EnsemblId, ids_dna, drop = FALSE]
    y <- as.numeric(t(as.matrix(assay(y_gene))))

    fold_r2 <- tapply(seq_along(foldid), foldid, function(idx) {
      cor(y[idx], preds[idx])^2
    })

    data.frame(
      fold = as.integer(names(fold_r2)),
      R2 = as.numeric(fold_r2),
      group = group,
      Symbol = symbol
    )
  }

  #per gene: paired Wilcoxon tests on fold-wise R2 between models, plus mean and sd of fold R2
  #genes that error are dropped (NULL); produces fold_results
  genes <- c(merged_df$Gene)

  results_list <- lapply(genes, function(symbol) {
    tryCatch({
      fold_DNA  <- extract_fold_r2(rpkm_se, "Genotype", symbol)
      fold_DNAm <- extract_fold_r2(rpkm_se, "DNAm", symbol)
      fold_meta <- extract_fold_r2(rpkm_se, "meta", symbol)
      fold_all <- bind_rows(fold_DNA, fold_DNAm, fold_meta)
      fold_wide <- fold_all %>% select(fold, group, R2) %>% pivot_wider(names_from = group, values_from = R2)
      w_meta_vs_DNA  <- wilcox.test(fold_wide$meta, fold_wide$Genotype, paired = TRUE, exact = TRUE)
      w_meta_vs_DNAm <- wilcox.test(fold_wide$meta, fold_wide$DNAm, paired = TRUE, exact = TRUE)
      w_DNA_vs_DNAm <- wilcox.test(fold_wide$Genotype, fold_wide$DNAm, paired = TRUE, exact = TRUE)
      data.frame(
        Gene = symbol,
        p_meta_vs_Genotype = w_meta_vs_DNA$p.value,
        p_meta_vs_DNAm     = w_meta_vs_DNAm$p.value,
        p_DNA_vs_DNAm     = w_DNA_vs_DNAm$p.value,
        mean_R2_meta  = mean(fold_wide$meta),
        mean_R2_DNA   = mean(fold_wide$Genotype),
        mean_R2_DNAm  = mean(fold_wide$DNAm),
        sd_R2_meta    = sd(fold_wide$meta),
        sd_R2_DNA     = sd(fold_wide$Genotype),
        sd_R2_DNAm    = sd(fold_wide$DNAm)
      )
    }, error = function(e) NULL)
  })
  warnings()
  fold_results <- bind_rows(results_list)

  #p value correction (FDR across all three comparisons together)
  all_p <- c(
    fold_results$p_meta_vs_Genotype,
    fold_results$p_meta_vs_DNAm,
    fold_results$p_DNA_vs_DNAm
  )

  all_p_fdr <- p.adjust(all_p, method = "fdr")

  #split back into columns
  n <- nrow(fold_results)
  fold_results$p_meta_vs_Genotype_FDR <- all_p_fdr[1:n]
  fold_results$p_meta_vs_DNAm_FDR      <- all_p_fdr[(n+1):(2*n)]
  fold_results$p_DNA_vs_DNAm_FDR       <- all_p_fdr[(2*n+1):(3*n)]

  #save output
  save(fold_results, file = file.path(base_dir, paste0(prefix, "all_results_gene_stats.RData")))

  #export as csv
  write.csv(fold_results, paste0(prefix, "all_results_gene_stats.csv"), row.names = FALSE)
}
