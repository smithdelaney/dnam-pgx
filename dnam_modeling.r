#title: DNAm_model
#author: Delaney Smith
#date: '2025-6-17'
#function: Estimate variance explained epigenetically in gene expression
#output: Files in Output 4 directory 'DNAm_model_output'

#load required packages
library(SummarizedExperiment)
library(tidyverse)
library(minfi)
library(data.table)
library(dplyr)
library(glmnet)
library(stringr)

#expression data: 'rpkm_pgkb' (object rpkm_se) for highly expressed genes, 'rpkm_pgx' for PGx genes
load(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/rpkm_pgkb.RData")
rpkm_normal <- rpkm_se
load(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/rpkm_pgx.RData")

#gene lists: 'filtered_genes' (highly expressed pharmGKB annotated genes) and 'pgx' (PGx high profile genes, regardless of expression)
load(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/filtered_genes.RData")
load(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/pgx_genes.RData")

#model both gene sets (a gene in both is modeled once, as PGx)
gene_cols <- intersect(colnames(pgx), colnames(filtered_genes))
filtered_genes <- rbind(cbind(pgx[, gene_cols], gene_set = "pgx"), cbind(filtered_genes[, gene_cols], gene_set = "normal"))
filtered_genes <- filtered_genes[!duplicated(filtered_genes$Symbol), ]

#remove x/y genes (applies to both sets)
filtered_genes <- filtered_genes[!tolower(filtered_genes$Chromosome) %in% c("chrx", "chry"), ]

#skip genes already modeled
remove_processed_genes <- function(filtered_genes, dir_path = ".") {
  processed_files <- list.files(dir_path, pattern = "\\.RData$", full.names = FALSE)
  keep <- !sapply(filtered_genes$Symbol, function(gene) {any(grepl(paste0("^", gene, "_DNAm\\.RData$"), processed_files))})
  filtered_genes[keep, ]
}

filtered_genes <- remove_processed_genes(filtered_genes, "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/DNAm_model_output")
#one gene per run (run_dnam_1x.sh calls this script repeatedly)
filtered_genes <- filtered_genes[1,]

#set lambda rule
lambda.rule <- "lambda.min"

#start loop through genes
for (i in 1:nrow(filtered_genes)) {

  #expression data for this gene's set
  rpkm_se <- if (filtered_genes$gene_set[i] == "pgx") rpkm_pgx else rpkm_normal

  #extract chromosome for gene
  chr <- str_extract(filtered_genes$Chromosome[i], "\\d+")
  chr_name <- paste0("chr", chr)

  #load DNAm data
  mfile_path <- file.path("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_1", paste0("mvalues_BIOS_", chr_name, ".RData"))
  load(mfile_path)

  #prepare covariates (DNAm covariates plus flowcell from mRNA)
  cvrt <- as.data.frame(colData(mvalues_chr))[, c("biobank_id", "sentrix_position","sample_plate", "sampling_age", "sex", "Mono_predicted", "Neu_predicted", "NK_predicted", "Eos_predicted", "CD8Tnv_predicted", "Treg_predicted", "Bnv_predicted", "Bmem_predicted" , "CD4Tmem_predicted", "CD4Tmem_predicted", "Baso_predicted","CD4Tnv_predicted"), drop = FALSE]
  e_cvrt <- as.data.frame(colData(rpkm_se))[, c("flowcell_num"), drop = FALSE]
  cvrt$biobank_id <- factor(cvrt$biobank_id)
  cvrt$sex <- factor(cvrt$sex)
  cvrt$sentrix_position<-factor(cvrt$sentrix_position)
  cvrt$sample_plate<-factor(cvrt$sample_plate)
  e_cvrt$flowcell_num<-factor(e_cvrt$flowcell_num)
  all_cvrt <- merge(cvrt, e_cvrt, by = 0)

  #subset all datasets to matching patients (in-loop since DNAm patients vary)
  common_columns <- Reduce(intersect, list(colnames(mvalues_chr), colnames(rpkm_se)))
  mvalues_chr<- mvalues_chr[, common_columns] #2524 samples
  rpkm_se <- rpkm_se[, common_columns] #2524 samples
  cvrt <- all_cvrt[all_cvrt$Row.names %in% common_columns, ] #2524 samples
  cvrt <- cvrt[, colnames(cvrt) != "Row.names"]
  cvrt<-model.matrix(~.,data=cvrt)[,-1]

  #subset the mRNA data by gene
  y <- rpkm_se[rownames(rpkm_se) == filtered_genes$Ensembl.Id[i], ]
  y <- t(as.matrix(assay(y)))

  #stratified 10 fold id generation based on sex and biobank id
  set.seed(1234)

  meta <- as.data.frame(colData(rpkm_se))
  meta$strata <- interaction(meta$sex, meta$biobank_id, drop = TRUE)
  foldid <- numeric(nrow(meta))

  for (s in unique(meta$strata)) {
    idx <- which(meta$strata == s)
    foldid[idx] <- sample(rep(1:10, length.out = length(idx)))
  }

  ######### DNAm ONLY MODEL ###############
  #format input matrix
  x <- t(as.matrix(assay(mvalues_chr)))
  x <- cbind(x, cvrt)

  cv.fit <- cv.glmnet(x, y, keep=TRUE, alpha=0.5, foldid=foldid)

  #extract cross-validated R2 at lambda.min
  id<-which(cv.fit$lambda == cv.fit$lambda.min)
  if(var(cv.fit$fit.preval[,id])!=0) {
    R2 <- cor(cv.fit$fit.preval[,id],y)^2
  } else {
    R2 <- 0
  }

  #save: cv model, R2 for CV, gene name
  DNAm_results <- list(EnsemblId = filtered_genes$Ensembl.Id[i], R2 = R2, cvfit = cv.fit)

  #free memory
  rm(cv.fit, x)
  gc()

  save(DNAm_results, file = paste0("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/DNAm_model_output/", filtered_genes$Symbol[i], "_DNAm.RData"))

}