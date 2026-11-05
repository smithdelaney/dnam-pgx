#title: analysis
#author: Delaney Smith
#date: '2025-6-25'
#function: summary statistics and figures for per-gene R2 of the DNA, DNAm and combined models
#output: printed summary stats and PNG figures in /home/dasmith/figures (pgx_ prefix for PGx genes, none for the rest)

#load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(tidyverse)
library(minfi)
library(data.table)
library(glmnet)
library(stringr)
library(ggrepel)

#directory of model summaries
base_dir <- "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4"

#run for the non-PGx genes (no file prefix) and the PGx genes (pgx_ prefix)
for (gene_set in c("normal", "pgx")) {

  prefix <- if (gene_set == "pgx") "pgx_" else ""
  cat("Gene set:", gene_set, "\n")

  #load results (merged_df), per-gene model comparison stats (fold_results, from post_process.r) and expression (y) data
  load(file = file.path(base_dir, paste0(prefix, "all_results_final.RData")))
  load(file = file.path(base_dir, paste0(prefix, "all_results_gene_stats.RData")))
  if (gene_set == "pgx") {
    load(file.path(base_dir, "rpkm_pgx.RData"))
    rpkm_se <- rpkm_pgx
  } else {
    load(file.path(base_dir, "rpkm_pgkb.RData"))
  }

  #a loop does not autoprint, so evaluate the analysis with autoprint as at top level
  withAutoprint({

    ############# per gene summary #############

    #results for a single gene
    merged_df %>%
      filter(Gene == "ABCA2")

    ############# top X genes #############

    #set top x genes
    X <- 5

    #top X by DNA model R2
    top_DNA <- merged_df %>%
      arrange(desc(R2_DNA)) %>%
      slice_head(n = X)
    print(top_DNA)

    #top X by DNAm model R2
    top_DNAm <- merged_df %>%
      arrange(desc(R2_DNAm)) %>%
      slice_head(n = X)
    print(top_DNAm)

    #top X by combined model R2
    top_comb <- merged_df %>%
      arrange(desc(R2_comb)) %>%
      slice_head(n = X)
    print(top_comb)

    ############# top X genes by maximum differences in R2 #############

    #rank by DNAm vs Genotype
    ranked_genes <- merged_df %>%
      mutate(DNAm_minus_Genotype = R2_DNAm - R2_DNA) %>%
      arrange(desc(DNAm_minus_Genotype))
    head(ranked_genes, X)

    #rank by combined vs best single (DNA or DNAm)
    ranked_genes_meta <- merged_df %>%
      rowwise() %>%
      mutate(
        best_single = max(R2_DNA, R2_DNAm, na.rm = TRUE),
        Meta_minus_best = R2_comb - best_single
      ) %>%
      ungroup() %>%
      arrange(desc(Meta_minus_best))
    head(ranked_genes_meta, X)

    ############## model level summary stats ##############

    #mean R2 and % of genes with R2 above 0.5, per model

    #DNA
    mean_R2 <- mean(merged_df$R2_DNA)
    n_above_0_5 <- sum(merged_df$R2_DNA > 0.5)
    percent <- (n_above_0_5 / nrow(merged_df)) * 100
    cat(sprintf("Mean R2 DNA = %.4f | %% > 0.5 = %.1f%%\n", mean_R2, percent))

    #DNAm
    mean_R2 <- mean(merged_df$R2_DNAm)
    n_above_0_5 <- sum(merged_df$R2_DNAm > 0.5)
    percent <- (n_above_0_5 / nrow(merged_df)) * 100
    cat(sprintf("Mean R2 DNAm = %.4f | %% > 0.5 = %.1f%%\n", mean_R2, percent))

    #comb
    mean_R2 <- mean(merged_df$R2_comb)
    n_above_0_5 <- sum(merged_df$R2_comb > 0.5)
    percent <- (n_above_0_5 / nrow(merged_df)) * 100
    cat(sprintf("Mean R2 comb = %.4f | %% > 0.5 = %.1f%%\n", mean_R2, percent))

    #Friedman test (non-parametric test for any difference between the three models)
    long_df <- merged_df %>%
      select(EnsemblId, R2_DNAm, R2_DNA, R2_comb) %>%
      pivot_longer(cols = starts_with("R2"), names_to = "Model", values_to = "R2")
    friedman.test(R2 ~ Model | EnsemblId, data = long_df)

    #follow up paired Wilcoxon tests for each model pair, BH adjusted
    models <- unique(long_df$Model)

    pairwise_wilcox <- combn(models, 2, function(m) {
      d1 <- long_df[long_df$Model == m[1], c("EnsemblId", "R2")]
      d2 <- long_df[long_df$Model == m[2], c("EnsemblId", "R2")]
      d  <- merge(d1, d2, by = "EnsemblId")

      data.frame(
        model1 = m[1],
        model2 = m[2],
        p.value = wilcox.test(d$R2.x, d$R2.y, paired = TRUE)$p.value
      )
    }, simplify = FALSE)

    pairwise_wilcox <- do.call(rbind, pairwise_wilcox)
    pairwise_wilcox$p.adj <- p.adjust(pairwise_wilcox$p.value, method = "BH")

    pairwise_wilcox %>%
      arrange(p.adj) %>%
      print(row.names = FALSE)

    #average gain in R2 between two modalities
    mean(merged_df$R2_comb - merged_df$R2_DNA, na.rm = TRUE)
    mean(merged_df$R2_comb - merged_df$R2_DNAm, na.rm = TRUE)
    mean(merged_df$R2_DNA - merged_df$R2_DNAm, na.rm = TRUE)

    #average proportion of combined model weight that is DNAm
    mean(merged_df$prop_dnam, na.rm = TRUE)

    ############## gene level summary stats ##############

    #fold-wise results for a specific gene
    fold_results %>%
      filter(Gene == "ADA")

    #number of genes with higher R2 for DNA, and how many of those are significant
    dna_better <- fold_results$mean_R2_DNA > fold_results$mean_R2_DNAm
    sum(dna_better)

    dna_better_sig <- dna_better & (fold_results$p_DNA_vs_DNAm_FDR < 0.05)
    sum(dna_better_sig)

    #and vice versa
    dnam_better <- fold_results$mean_R2_DNAm > fold_results$mean_R2_DNA
    sum(dnam_better)

    dnam_better_sig <- dnam_better & (fold_results$p_DNA_vs_DNAm_FDR < 0.05)
    sum(dnam_better_sig)

    #comb vs best single modality
    best_R2 <- pmax(fold_results$mean_R2_DNA, fold_results$mean_R2_DNAm)

    meta_beats_best <- fold_results$mean_R2_meta > best_R2

    meta_beats_best_sig <- meta_beats_best & (
      (fold_results$mean_R2_DNA >= fold_results$mean_R2_DNAm & fold_results$p_meta_vs_Genotype_FDR < 0.05) |
      (fold_results$mean_R2_DNAm >  fold_results$mean_R2_DNA  & fold_results$p_meta_vs_DNAm_FDR     < 0.05)
    )

    sum(meta_beats_best)
    sum(meta_beats_best_sig)

    #comb vs DNA
    comb_better <- fold_results$mean_R2_meta > fold_results$mean_R2_DNA
    sum(comb_better)

    comb_better_sig <- comb_better & (fold_results$p_meta_vs_Genotype_FDR < 0.05)
    sum(comb_better_sig)

    ###################### HLA and top x gene sets ######################

    #gene sets (hard-coded from results)
    hla_genes <- c("HLA-B", "HLA-A", "HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1", "HLA-DRB5")
    dna_genes <- c("ERAP2", "IGHMBP2", "SH3BP2", "CAPN10", "EGLN2")
    dnam_genes <- c("HLA-DQB1", "TBX21", "APBB1", "HLA-DQA1", "TCL1A")
    comb_genes <- c("SULT1A1", "SLC29A1", "IL6R", "FAAH", "ADA")

    #mean R2 of each model in HLA genes
    mean(merged_df$R2_DNAm[merged_df$Gene %in% hla_genes], na.rm = TRUE)
    mean(merged_df$R2_DNA[merged_df$Gene %in% hla_genes], na.rm = TRUE)
    mean(merged_df$R2_comb[merged_df$Gene %in% hla_genes], na.rm = TRUE)

    fold_results %>%
      filter(Gene %in% hla_genes)

    #mean R2 and DNAm weight in the top DNA genes
    c(
      avg_R2_DNA    = mean(merged_df$R2_DNA[merged_df$Gene %in% dna_genes], na.rm = TRUE),
      avg_prop_DNAm = mean(merged_df$prop_dnam[merged_df$Gene %in% dna_genes], na.rm = TRUE)
    )

    #mean R2 and DNAm weight in the top DNAm genes
    c(
      avg_R2_DNAm    = mean(merged_df$R2_DNAm[merged_df$Gene %in% dnam_genes], na.rm = TRUE),
      avg_prop_DNAm = mean(merged_df$prop_dnam[merged_df$Gene %in% dnam_genes], na.rm = TRUE)
    )

    #mean DNAm weight in HLA and top combined genes
    print(mean(merged_df$prop_dnam[merged_df$Gene %in% hla_genes], na.rm = TRUE))
    print(mean(merged_df$prop_dnam[merged_df$Gene %in% comb_genes], na.rm = TRUE))

    ############# plotting #############

    #output directory for plots
    plot_dir <- "/home/dasmith/figures"

    #set theme
    my_theme <- theme_minimal(base_size = 16) +
      theme(
        axis.text = element_text(size = 16, face = "plain"),
        axis.title = element_text(size = 16, face = "plain"),
        plot.title = element_text(hjust = 0.5, face = "plain"),
        panel.grid.minor = element_blank(),
      )

    theme_set(my_theme)

    #color refs
    #"#89CFF0",  # pastel blue, DNA
    #"#77DD77",      # pastel green, DNAm
    #"#FF6961"       # pastel red, comb

    #R2 histogram for one model (adjust model and color as needed)
    model <- "DNA"  # or "DNAm", "comb"

    p <- ggplot(merged_df, aes(x = .data[[paste0("R2_", model)]])) +
      geom_histogram(bins = 30, fill = "#89CFF0", color = "white") +
      labs(
        title = paste(model),
        x = "R²",
        y = "Count"
      )

    fig_path <- file.path(plot_dir, paste0(prefix, "final_R2_histogram_", model, ".png"))
    ggsave(fig_path, plot = p, width = 4, height = 3, dpi = 300)

    #overlapping R2 density plots

    #long format, one row per gene and model
    results_df <- merged_df %>%
      select(Gene, R2_DNA, R2_DNAm, R2_comb) %>%
      pivot_longer(
        cols = starts_with("R2_"),
        names_to = "group",
        values_to = "R2"
      ) %>%
      mutate(
        group = recode(
          group,
          R2_DNA  = "Genotype",
          R2_DNAm = "DNAm",
          R2_comb = "meta"
        )
      )

    group_colors <- c(
      Genotype = "#89CFF0",
      DNAm     = "#77DD77",
      meta     = "#FF6961"
    )

    #DNA
    p_dna <- ggplot(
      results_df %>% filter(group == "Genotype"),
      aes(x = R2, fill = group)
    ) +
      geom_density(alpha = 0.6) +
      scale_fill_manual(values = group_colors) +
      scale_x_continuous(expand = expansion(mult = 0)) +
      scale_y_continuous(expand = expansion(mult = 0)) +
      labs(title = "", x = "R²", y = "Density")

    ggsave(
      file.path(plot_dir, paste0(prefix, "final_R2_density_DNA.png")),
      plot = p_dna, width = 8, height = 6, dpi = 300
    )

    #DNA + DNAm
    p_dna_dnam <- ggplot(
      results_df %>% filter(group %in% c("Genotype", "DNAm")),
      aes(x = R2, fill = group)
    ) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = group_colors) +
      scale_x_continuous(expand = expansion(mult = 0)) +
      scale_y_continuous(expand = expansion(mult = 0)) +
      labs(title = "", x = "R²", y = "Density", fill = "Group")

    ggsave(
      file.path(plot_dir, paste0(prefix, "final_R2_density_DNA_DNAm.png")),
      plot = p_dna_dnam, width = 8, height = 6, dpi = 300
    )

    #all three
    p_all <- ggplot(
      results_df,
      aes(x = R2, fill = group)
    ) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = group_colors) +
      scale_x_continuous(expand = expansion(mult = 0)) +
      scale_y_continuous(expand = expansion(mult = 0)) +
      labs(title = "", x = "R²", y = "Density", fill = "Group")

    ggsave(
      file.path(plot_dir, paste0(prefix, "final_R2_density_all.png")),
      plot = p_all, width = 8, height = 8, dpi = 300
    )

    #R2 correlation scatter plots, plain and gene-labeled versions

    #DNAm vs DNA
    #lm_fit <- lm(R2_DNA ~ R2_DNAm, data = merged_df)
    #lm_r2 <- summary(lm_fit)$r.squared
    #print(round(lm_r2, 3))

    spearman_corr <- cor(merged_df$R2_DNAm, merged_df$R2_DNA, method = "spearman", use = "complete.obs")
    print(round(spearman_corr, 3))

    #plot
    d <- ggplot(merged_df, aes(x = R2_DNAm, y = R2_DNA)) +
      geom_point(color = "#7f7f7f", size = 4, alpha = 0.7) +  #larger, colored dots
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
      geom_smooth(method = "lm", se = FALSE, color = "black") +
      labs(
        title = "",
        x = "R² DNAm",
        y = "R² DNA"
      )
    #save
    ggsave(file.path(plot_dir, paste0(prefix, "final_DNAm_DNA_corr.png")), plot = d, width = 8, height = 6, dpi = 300)

    #make labeled version
    d_labels <- d +
      geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 50)  #adjust size and max.overlaps as needed
    ggsave(file.path(plot_dir, paste0(prefix, "DNAm_DNA_corr_labeled.png")), plot = d_labels, width = 10, height = 8, dpi = 300)

    #DNAm vs combined
    spearman_corr <- cor(merged_df$R2_DNAm, merged_df$R2_comb, method = "spearman", use = "complete.obs")
    print(round(spearman_corr, 3))

    #plot
    d <- ggplot(merged_df, aes(x = R2_DNAm, y = R2_comb)) +
      geom_point(color = "#33a02c", size = 4, alpha = 0.7) +  #larger, colored dots
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
      geom_smooth(method = "lm", se = FALSE, color = "black") +
      labs(
        title = "",
        x = "R² DNAm",
        y = "R² Combined"
      )
    #save
    ggsave(file.path(plot_dir, paste0(prefix, "final_DNAm_comb_corr.png")), plot = d, width = 8, height = 6, dpi = 300)

    #make labeled version
    d_labels <- d +
      geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 50)  #adjust size and max.overlaps as needed
    ggsave(file.path(plot_dir, paste0(prefix, "DNAm_comb_corr_labeled.png")), plot = d_labels, width = 10, height = 8, dpi = 300)

    #DNA vs combined
    spearman_corr <- cor(merged_df$R2_DNA, merged_df$R2_comb, method = "spearman", use = "complete.obs")
    print(round(spearman_corr, 3))

    #plot
    d <- ggplot(merged_df, aes(x = R2_DNA, y = R2_comb)) +
      geom_point(color = "#1f78b4", size = 4, alpha = 0.7) +  #larger, colored dots
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
      geom_smooth(method = "lm", se = FALSE, color = "black") +
      labs(
        title = "",
        x = "R² DNA",
        y = "R² Combined"
      )
    #save
    ggsave(file.path(plot_dir, paste0(prefix, "final_DNA_comb_corr.png")), plot = d, width = 8, height = 6, dpi = 300)

    #make labeled version
    d_labels <- d +
      geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 50)  #adjust size and max.overlaps as needed
    ggsave(file.path(plot_dir, paste0(prefix, "DNA_comb_corr_labeled.png")), plot = d_labels, width = 10, height = 8, dpi = 300)

    #overall summary figure: mean R2 per model with 95% CI

    #long format, one row per gene and model
    long_df <- merged_df %>%
      select(R2_DNA, R2_DNAm, R2_comb) %>%
      pivot_longer(cols = everything(), names_to = "Group", values_to = "R2")

    summary_stats <- long_df %>%
      group_by(Group) %>%
      summarise(
        mean_R2  = mean(R2, na.rm = TRUE),
        sd_R2    = sd(R2, na.rm = TRUE),
        n        = sum(!is.na(R2)),
        se       = sd(R2, na.rm = TRUE) / sqrt(sum(!is.na(R2))),
        ci_lower = mean(R2, na.rm = TRUE) - 1.96 * se,
        ci_upper = mean(R2, na.rm = TRUE) + 1.96 * se,
        .groups = "drop"
      )

    #custom bar labels
    bar_labels <- c(
      "R2_DNA"  = "DNA",
      "R2_DNAm" = "DNAm",
      "R2_comb" = "Combined"
    )

    #set order
    summary_stats$Group <- factor(
      summary_stats$Group,
      levels = c("R2_DNAm", "R2_DNA", "R2_comb")
    )

    #plot
    s <- ggplot(summary_stats, aes(x = Group, y = mean_R2, fill = Group)) +
      geom_bar(stat = "identity", color = "black", width = 0.6) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
      scale_fill_manual(values = c(
        "R2_DNA" = "#89CFF0",    #pastel blue
        "R2_DNAm" = "#77DD77",   #pastel green
        "R2_comb" = "#FF6961"    #pastel red
      )) +
      theme(
        panel.grid.major.x = element_blank(),  #remove vertical major grid lines
        panel.grid.minor.x = element_blank()   #remove vertical minor grid lines
      )+
      scale_x_discrete(labels = bar_labels) +  #set custom bar names
      ylab("R²") +
      xlab("")

    ggsave(file.path(plot_dir, paste0(prefix, "final_summary_bars.png")), plot = s, width = 8, height = 6, dpi = 300)

    ###################### PGx-specific stats ######################

    if (gene_set == "pgx") {

      #average RPKM for all genes except HLA-B (included in other set)
      print(mean(assay(rpkm_se)[rownames(rpkm_se) != "ENSG00000234745", ]))

      #average RPKM for 'null' genes
      genes_of_interest <- c(
        "ENSG00000197408",  #CYP2B6
        "ENSG00000138109",  #CYP2C9
        "ENSG00000165841",  #CYP2C19
        "ENSG00000134538",  #SLCO1B1
        "ENSG00000241635"   #UGT1A1
      )
      print(mean(assay(rpkm_se)[genes_of_interest, , drop = FALSE]))
    }

  }, echo = FALSE)
}
