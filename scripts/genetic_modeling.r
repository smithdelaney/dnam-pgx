#title: genetic_modelling_1x
#author: Delaney Smith
#date: '2025-6-17'
#function: Estimate variance explained genetically in gene expression
#output: Files in Output 4 directory

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

#object is called 'counts' --> raw mRNA data w/covariate info
load(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_3/Output_03_mRNA_VKORC1_c.RData")

#remove x/y chromosome genes (applies to both sets)
filtered_genes <- filtered_genes[!tolower(filtered_genes$Chromosome) %in% c("chrx", "chry"), ]

#skip genes already modeled
remove_processed_genes <- function(filtered_genes, dir_path = ".") {
  processed_files <- list.files(dir_path, pattern = "\\.RData$", full.names = FALSE)
  keep <- !sapply(filtered_genes$Symbol, function(gene) {any(grepl(paste0("^", gene, "_DNA\\.RData$"), processed_files))})
  filtered_genes[keep, ]
}

filtered_genes <- remove_processed_genes(filtered_genes, "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/DNA_model_output")
#one gene per run for memory (run_DNA_1x.sh calls this script repeatedly)
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

  #load DNA data, object called chr{chr}_BIOS_dosages
  file_path <- file.path("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2", paste0(chr_name, "_BIOS_dosages.RData"))
  load(file_path)
  name <- paste0(chr_name, "_BIOS_dosages")
  dosages <- get(name)

  #prepare covariates
  cvrt <- as.data.frame(colData(counts))[, c("biobank_id", "flowcell_num", "sampling_age", "sex", "Mono_predicted", "Neu_predicted", "NK_predicted", "Eos_predicted", "CD8Tnv_predicted", "Treg_predicted", "Bnv_predicted", "Bmem_predicted" , "CD4Tmem_predicted", "CD4Tmem_predicted", "Baso_predicted","CD4Tnv_predicted"), drop = FALSE]
  cvrt$biobank_id <- factor(cvrt$biobank_id)
  cvrt$sex <- factor(cvrt$sex)
  cvrt$flowcell_num<-factor(cvrt$flowcell_num)

  #subset all datasets to matching patients (in-loop since DNA patients vary)
  common_columns <- Reduce(intersect, list(colnames(dosages), colnames(rpkm_se)))
  dosages <- dosages[, common_columns] #2524 samples
  rpkm_se <- rpkm_se[, common_columns] #2524 samples
  #cvrt <- cvrt[cvrt$Row.names %in% common_columns, ] #2524 samples
  #cvrt <- cvrt[, colnames(cvrt) != "Row.names"]
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

  ######### DNA ONLY MODEL ###############

  #format input matrix (dosages plus covariates)
  x <- t(as.matrix(dosages))
  x <- cbind(x, cvrt)

  #fit cross-validated elastic net (DNA only)
  cv.fit <- cv.glmnet(x, y, keep=TRUE, alpha=0.5, foldid=foldid)

  #extract cross-validated R2 at lambda.min
  id<-which(cv.fit$lambda == cv.fit$lambda.min)
  if(var(cv.fit$fit.preval[,id])!=0) {
    R2 <- cor(cv.fit$fit.preval[,id],y)^2
  } else {
    R2 <- 0
  }

  #save: cv model, R2 for CV, gene name
  DNA_results <- list(EnsemblId = filtered_genes$Ensembl.Id[i],R2 = R2,cvfit = cv.fit)

  #free memory
  rm(cv.fit, x)
  gc()

  save(DNA_results, file = paste0("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/DNA_model_output/", filtered_genes$Symbol[i], "_DNA.RData"))

}
