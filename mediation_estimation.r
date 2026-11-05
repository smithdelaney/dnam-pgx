#title: Mediation
#author: Delaney Smith
#date: '2026-1-5'
#function: mediation analysis of how much DNA model performance is mediated by DNAm, using out-of-fold predictions with bootstrap CIs
#output: mediation.RData and pgx_mediation.RData (mediation_df for non-PGx and PGx genes) in Output 4 directory

#load libraries (check for not needed packages)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(tidyverse)
library(minfi)
library(data.table)
library(glmnet)
library(stringr)

#directories
base_dir <- "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4"
dna_dir  <- file.path(base_dir, "DNA_model_output")
dnam_dir <- file.path(base_dir, "DNAm_model_output")

######### gene sets #################

hla_genes <- c("HLA-B", "HLA-A", "HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1", "HLA-DRB5")
dna_genes <- c("ERAP2", "IGHMBP2", "SH3BP2", "CAPN10", "EGLN2")
dnam_genes <- c("HLA-DQB1", "TBX21", "APBB1", "HLA-DQA1", "TCL1A")
comb_genes <- c("SULT1A1", "SLC29A1", "IL6R", "FAAH", "ADA")

############# mediation analysis #############

#mediation function: DNAm as mediator of the DNA effect on expression
#uses out-of-fold predictions of the DNA and DNAm cv.glmnet fits
#returns list of R2 (total, direct, mediated), proportion mediated, coefficient effects and optional bootstrap draws/CIs
mediation <- function(
  dna_cvfit,
  dnam_cvfit,
  y,
  ids = NULL,
  lambda = c("lambda.min", "lambda.1se"),
  bootstrap = FALSE,
  n_boot = 500,
  seed = 1
) {
  lambda <- match.arg(lambda)

  #extract out-of-fold preds
  get_oof_pred <- function(cvfit, lambda_choice) {
    if (is.null(cvfit$fit.preval)) stop("cvfit$fit.preval is NULL. Refit with keep=TRUE.")
    lam <- cvfit[[lambda_choice]]
    j <- which.min(abs(cvfit$glmnet.fit$lambda - lam))
    pred <- as.numeric(cvfit$fit.preval[, j])
    sid  <- rownames(cvfit$fit.preval)
    if (is.null(sid)) stop("rownames(cvfit$fit.preval) are NULL; need sample IDs to align.")
    list(pred = pred, ids = sid)
  }

  #align y to sample IDs
  if (is.null(ids)) {
    if (is.null(names(y))) stop("Provide `ids` or set names(y) to sample IDs.")
    ids_y <- names(y)
  } else {
    ids_y <- ids
  }
  if (length(y) != length(ids_y)) stop("length(y) must equal length(ids).")

  #pull predictions
  p_dna  <- get_oof_pred(dna_cvfit,  lambda)
  p_dnam <- get_oof_pred(dnam_cvfit, lambda)

  common_ids <- Reduce(intersect, list(p_dna$ids, p_dnam$ids, ids_y))
  if (length(common_ids) < 10) stop("Too few overlapping samples after alignment.")

  DNA  <- p_dna$pred[match(common_ids, p_dna$ids)]
  DNAm <- p_dnam$pred[match(common_ids, p_dnam$ids)]
  Y <- as.numeric(y[match(common_ids, ids_y)])

  df <- data.frame(Y = Y, DNA = as.numeric(DNA), DNAm = as.numeric(DNAm))

  #mediation set up

  #mediator
  fit_M <- lm(DNAm ~ DNA, data = df)
  a1 <- unname(coef(fit_M)["DNA"])

  #DNA only
  fit_A <- lm(Y ~ DNA, data = df)
  R2_DNA_total <- summary(fit_A)$r.squared

  #DNA and DNAm (full model)
  fit_B <- lm(Y ~ DNA + DNAm, data = df)
  R2_B <- summary(fit_B)$r.squared

  #DNAm only (partial model)
  fit_B_noDNA <- lm(Y ~ DNAm, data = df)
  R2_B_noDNA <- summary(fit_B_noDNA)$r.squared

  #R2 DNA direct
  R2_DNA_direct <- R2_B - R2_B_noDNA

  #R2 mediated
  R2_DNA_mediated <- R2_DNA_total - R2_DNA_direct
  prop_DNA_mediated <- if (R2_DNA_total > 0) R2_DNA_mediated / R2_DNA_total else NA_real_

  #coefficient-scale effects
  b2 <- unname(coef(fit_B)["DNAm"])
  direct_coef <- unname(coef(fit_B)["DNA"])
  indirect_coef <- a1 * b2
  total_coef <- direct_coef + indirect_coef

  res <- list(
    n = length(common_ids),
    ids = common_ids,
    lambda = lambda,
    a1 = a1,
    b2 = b2,
    indirect_coef = indirect_coef,
    direct_coef = direct_coef,
    total_coef = total_coef,
    R2_DNA_total = R2_DNA_total,
    R2_DNA_direct = R2_DNA_direct,
    R2_DNA_mediated = R2_DNA_mediated,
    prop_DNA_mediated = prop_DNA_mediated
  )

  #bootstrapping (resample samples, refit the lm models, 2.5/50/97.5% quantiles)
  if (isTRUE(bootstrap)) {
    set.seed(seed)
    boot <- matrix(NA_real_, nrow = n_boot, ncol = 5)
    colnames(boot) <- c(
      "indirect_coef", "direct_coef",
      "R2_DNA_direct", "R2_DNA_mediated", "prop_DNA_mediated"
    )

    for (b in seq_len(n_boot)) {
      idx <- sample.int(nrow(df), replace = TRUE)
      dfi <- df[idx, , drop = FALSE]

      #mediator model: DNAm ~ DNA
      fit_Mb <- lm(DNAm ~ DNA, data = dfi)
      a1b <- unname(coef(fit_Mb)["DNA"])

      #total DNA signal
      fit_Ab <- lm(Y ~ DNA, data = dfi)
      R2_Ab <- summary(fit_Ab)$r.squared

      #outcome model with mediator
      fit_Bb <- lm(Y ~ DNA + DNAm, data = dfi)
      fit_Bb_noDNA <- lm(Y ~ DNAm, data = dfi)

      #partial R2 of DNA given DNAm
      R2_direct_b <- summary(fit_Bb)$r.squared -
        summary(fit_Bb_noDNA)$r.squared

      #mediated DNA variance (CLAMPED)
      R2_med_b <- max(0, R2_Ab - R2_direct_b)
      prop_b <- if (R2_Ab > 0) R2_med_b / R2_Ab else NA_real_

      #coefficient-scale effects
      b2b <- unname(coef(fit_Bb)["DNAm"])
      direct_b <- unname(coef(fit_Bb)["DNA"])
      indirect_b <- a1b * b2b

      boot[b, ] <- c(indirect_b, direct_b, R2_direct_b, R2_med_b, prop_b)
    }

    ci <- apply(
      boot, 2, quantile,
      probs = c(0.025, 0.5, 0.975),
      na.rm = TRUE
    )

    res$bootstrap <- list(
      n_boot = n_boot,
      seed = seed,
      draws = boot,
      ci = ci
    )
  }

  res
}

#run for the non-PGx genes (no file prefix) and the PGx genes (pgx_ prefix)
mediation_results <- list()
for (gene_set in c("normal", "pgx")) {

  prefix <- if (gene_set == "pgx") "pgx_" else ""
  cat("Gene set:", gene_set, "\n")

  #load results (merged_df defines the gene set) and expression (y) data
  load(file = file.path(base_dir, paste0(prefix, "all_results_final.RData")))
  set_genes <- merged_df$Gene
  if (gene_set == "pgx") {
    load(file.path(base_dir, "rpkm_pgx.RData"))
    rpkm_se <- rpkm_pgx
  } else {
    load(file.path(base_dir, "rpkm_pgkb.RData"))
  }

  #get file paths
  dna_files <- list.files(dna_dir, pattern = "_DNA\\.RData$", full.names = TRUE)

  #per gene: load DNA and DNAm fits, run mediation with bootstrap
  #produces rows, a list of one-row data frames per gene
  rows <- list()
  for (dna_file in dna_files) {

    gene_name <- sub("_DNA\\.RData$", "", basename(dna_file))
    if (!gene_name %in% set_genes) next
    print(gene_name)

    dnam_file <- file.path(dnam_dir, paste0(gene_name, "_DNAm.RData"))
    if (!file.exists(dnam_file)) next

    #load DNA
    loaded_dna <- load(dna_file)
    if (!"DNA_results" %in% loaded_dna) {
      warning(
        "Skipping ", gene_name, ": DNA file does not contain DNA_results (has: ",
        paste(loaded_dna, collapse = ", "), ")"
      )
      next
    }
    d_cvfit <- DNA_results$cvfit

    #load DNAm
    loaded_dnam <- load(dnam_file)
    if (!"DNAm_results" %in% loaded_dnam) {
      warning(
        "Skipping ", gene_name, ": DNAm file does not contain DNAm_results (has: ",
        paste(loaded_dnam, collapse = ", "), ")"
      )
      next
    }
    m_cvfit <- DNAm_results$cvfit

    #expression vector
    y_gene <- rpkm_se[rownames(rpkm_se) == DNAm_results$EnsemblId, , drop = FALSE]
    if (nrow(y_gene) != 1) next

    y <- as.numeric(t(as.matrix(assay(y_gene))))
    names(y) <- colnames(assay(y_gene))

    #run mediation
    res <- mediation(
      dna_cvfit = d_cvfit,
      dnam_cvfit = m_cvfit,
      y = y,
      lambda = "lambda.min",
      bootstrap = TRUE,
      n_boot = 100
    )

    #store
    ci <- if (!is.null(res$bootstrap$ci)) res$bootstrap$ci else NULL
    rows[[gene_name]] <- data.frame(
      gene = gene_name,
      EnsemblId = DNAm_results$EnsemblId,
      n = res$n,
      R2_DNA_total = res$R2_DNA_total,
      R2_DNA_direct = res$R2_DNA_direct,
      R2_DNA_mediated = res$R2_DNA_mediated,
      prop_DNA_mediated = res$prop_DNA_mediated,
      indirect_coef = res$indirect_coef,
      direct_coef = res$direct_coef,
      #CIs (2.5%, 50%, 97.5%) for prop and mediated R2 from bootstrapping
      prop_ci_low  = if (!is.null(ci)) ci["2.5%",  "prop_DNA_mediated"] else NA_real_,
      prop_ci_med  = if (!is.null(ci)) ci["50%",   "prop_DNA_mediated"] else NA_real_,
      prop_ci_high = if (!is.null(ci)) ci["97.5%", "prop_DNA_mediated"] else NA_real_,
      R2med_ci_low  = if (!is.null(ci)) ci["2.5%",  "R2_DNA_mediated"] else NA_real_,
      R2med_ci_med  = if (!is.null(ci)) ci["50%",   "R2_DNA_mediated"] else NA_real_,
      R2med_ci_high = if (!is.null(ci)) ci["97.5%", "R2_DNA_mediated"] else NA_real_,
      stringsAsFactors = FALSE
    )

    rm(DNA_results, DNAm_results, d_cvfit, m_cvfit, y_gene, y, res)
    gc()
  }

  #combine per-gene rows into mediation_df
  mediation_df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

  #save
  save(mediation_df, file = file.path(base_dir, paste0(prefix, "mediation.RData")))

  mediation_results[[gene_set]] <- mediation_df
}

################ exploration ################

#gene lists below are non-PGx genes
mediation_df <- mediation_results$normal

#mean median proportion mediated across genes
mean(mediation_df$prop_ci_med)

#gene-wise
mediation_df$prop_ci_med[mediation_df$gene == "IGHMBP2"]

#top 5 dna
mean(mediation_df$prop_ci_med[mediation_df$gene %in% dna_genes])
range((mediation_df$prop_ci_med[mediation_df$gene %in% dna_genes]))

#top 5 dnam
mean(mediation_df$prop_ci_med[mediation_df$gene %in% dnam_genes])
range((mediation_df$prop_ci_med[mediation_df$gene %in% dnam_genes]))

#top 5 hla
mean(mediation_df$prop_ci_med[mediation_df$gene %in% hla_genes])
range((mediation_df$prop_ci_med[mediation_df$gene %in% hla_genes]))

#top 5 combined
mean(mediation_df$prop_ci_med[mediation_df$gene %in% comb_genes])
range((mediation_df$prop_ci_med[mediation_df$gene %in% comb_genes]))

#top 5 for mediation
mediation_df %>%
  arrange(desc(prop_ci_med)) %>%
  slice_head(n = 5)

#bottom 5 mediation
mediation_df %>%
  arrange(prop_ci_med) %>%
  slice_head(n = 5)
