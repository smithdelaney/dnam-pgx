#title: DNAm_mRNA_processing
#author: Delaney Smith, adapted from Yunfeng
#date: '2025-2-26'
#function: Prepare blood DNAm (mvalues, betas) and mRNA (RPKM) data for modeling
#output: Files in Output 1, 2 and 4 directories

#load packages
library(SummarizedExperiment)
library(FDb.InfiniumMethylation.hg19)
library(tidyverse)
library(DNAmArray)
library(minfi)
library(data.table)

#load DNA-methylation data (large, takes about 10 minutes)
load(file="/exports/archive/kft-pgx-archive/BIOS/methData_Mvalues_BIOS_Freeze2_unrelated.RData")
dim(mvalues)   # 481388    4386

#remove X and Y chromosome
mvalues<-dropSeqlevels(mvalues,c("chrX","chrY"),pruning.mode = "coarse")
dim(mvalues)   # 470211   4386

#load maskProbes data
maskProbes_450K <- as.data.frame(fread(file = "/home/dasmith/Input_01_HM450.hg19.manifest.tsv"))

#remove unreliable probes (60466 hg19 probes recommended for removal, Zhou et al.)
maskProbes_450K <- filter(maskProbes_450K,MASK_general=="TRUE")
mvalues_BIOS <- mvalues[!(rownames(mvalues) %in% maskProbes_450K$probeID),]
dim(mvalues_BIOS) # 412373   4386

#drop samples on sample_plate levels OV0192DNA001 and OV0192DNA002, which cause singularities in the models (Yunfeng et al.)
idx <- which(mvalues_BIOS$sample_plate == "OV0192DNA001" | mvalues_BIOS$sample_plate == "OV0192DNA002")
mvalues_BIOS <- mvalues_BIOS[,-idx]
dim(mvalues_BIOS)  # 412373   4045

#remove samples (N=263) with no age
table(is.na(mvalues_BIOS$sampling_age))
col.mvalues.BIOS<-as.data.frame(colData(mvalues_BIOS))
col.mvalues.BIOS<-col.mvalues.BIOS[!is.na(col.mvalues.BIOS$sampling_age),]
mvalues_BIOS<-mvalues_BIOS[,colnames(mvalues_BIOS) %in% rownames(col.mvalues.BIOS)]
dim(mvalues_BIOS)  # 412373   4123

#load beta values
load(file="/exports/archive/kft-pgx-archive/BIOS/methData_Betas_BIOS_Freeze2_unrelated.RData")
dim(betas)      # 481388   4386

#only keep beta values for probes in 450k annotation and already passing QC from the mvalues work above
betas_BIOS <- betas[!(rownames(betas) %in% maskProbes_450K$probeID),]
betas_BIOS<-betas_BIOS[,colnames(mvalues_BIOS)]
dim(betas_BIOS) # 422156   4045

#remove sex-mismatched or missing samples
predictedSex <- as.data.frame(getSex.DNAmArray(assay(betas_BIOS)))
colnames(predictedSex)<-"predictedSex"
predictedSex$predictedSex[grepl("Male",predictedSex$predictedSex)]<-"male"
predictedSex$predictedSex[grepl("Female",predictedSex$predictedSex)]<-"female"
mvalues_BIOS$predictedSex<-predictedSex$predictedSex
mvalues_BIOS$correctSex<-ifelse(colData(mvalues_BIOS)$predictedSex==colData(mvalues_BIOS)$sex,TRUE,FALSE)
table(colData(mvalues_BIOS)$correctSex) # FALSE:11  TRUE:4032 NA:2(sex is missing)

col.mvalues.BIOS<-as.data.frame(colData(mvalues_BIOS))
col.mvalues.BIOS<-filter(col.mvalues.BIOS,correctSex==TRUE)
mvalues_BIOS<-mvalues_BIOS[,rownames(col.mvalues.BIOS)]
betas_BIOS<-betas_BIOS[,rownames(col.mvalues.BIOS)]
betas_BIOS<-betas_BIOS[rownames(mvalues_BIOS),]
dim(mvalues_BIOS) # 412373  4032
dim(betas_BIOS)   # 412373  4032

#one problematic female sample ("BIOS691E2BA5") should be removed based on PCA from Li et al 2024.
mvalues_BIOS<-mvalues_BIOS[,!(colnames(mvalues_BIOS)%in%"BIOS691E2BA5")]
betas_BIOS<-betas_BIOS[,!(colnames(betas_BIOS)%in%"BIOS691E2BA5")]
dim(mvalues_BIOS) # 412373 4031
dim(betas_BIOS) # 412373 4031

#predicted cell counts for all BIOS data by IDOL-ext (Li et al.)
load(file="/home/dasmith/Input_01_EpiDISH_Cellcounts_BIOS.RData") #move file to storage
Cellcounts_BIOS_EpiDISH<-Cellcounts_BIOS_EpiDISH[colnames(mvalues_BIOS),]

#add predicted cell counts to betas and mvalues
betas_BIOS$CD4Tnv_predicted <- Cellcounts_BIOS_EpiDISH$CD4Tnv
betas_BIOS$Baso_predicted <- Cellcounts_BIOS_EpiDISH$Baso
betas_BIOS$CD4Tmem_predicted <- Cellcounts_BIOS_EpiDISH$CD4Tmem
betas_BIOS$Bmem_predicted <- Cellcounts_BIOS_EpiDISH$Bmem
betas_BIOS$Bnv_predicted <- Cellcounts_BIOS_EpiDISH$Bnv
betas_BIOS$Treg_predicted <- Cellcounts_BIOS_EpiDISH$Treg
betas_BIOS$CD8Tmem_predicted <- Cellcounts_BIOS_EpiDISH$CD8Tmem
betas_BIOS$CD8Tnv_predicted <- Cellcounts_BIOS_EpiDISH$CD8Tnv
betas_BIOS$Eos_predicted <- Cellcounts_BIOS_EpiDISH$Eos
betas_BIOS$NK_predicted <- Cellcounts_BIOS_EpiDISH$NK
betas_BIOS$Neu_predicted <- Cellcounts_BIOS_EpiDISH$Neu
betas_BIOS$Mono_predicted <- Cellcounts_BIOS_EpiDISH$Mono

mvalues_BIOS$CD4Tnv_predicted <- Cellcounts_BIOS_EpiDISH$CD4Tnv
mvalues_BIOS$Baso_predicted <- Cellcounts_BIOS_EpiDISH$Baso
mvalues_BIOS$CD4Tmem_predicted <- Cellcounts_BIOS_EpiDISH$CD4Tmem
mvalues_BIOS$Bmem_predicted <- Cellcounts_BIOS_EpiDISH$Bmem
mvalues_BIOS$Bnv_predicted <- Cellcounts_BIOS_EpiDISH$Bnv
mvalues_BIOS$Treg_predicted <- Cellcounts_BIOS_EpiDISH$Treg
mvalues_BIOS$CD8Tmem_predicted <- Cellcounts_BIOS_EpiDISH$CD8Tmem
mvalues_BIOS$CD8Tnv_predicted <- Cellcounts_BIOS_EpiDISH$CD8Tnv
mvalues_BIOS$Eos_predicted <- Cellcounts_BIOS_EpiDISH$Eos
mvalues_BIOS$NK_predicted <- Cellcounts_BIOS_EpiDISH$NK
mvalues_BIOS$Neu_predicted <- Cellcounts_BIOS_EpiDISH$Neu
mvalues_BIOS$Mono_predicted <- Cellcounts_BIOS_EpiDISH$Mono

#save whole files (mvalues to Output_2, betas to the Imputation root)
save(mvalues_BIOS,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/methData_mvalues_BIOS.RData")
save(betas_BIOS,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/methData_betas_BIOS.RData")

#reload mvalues_BIOS (from Output_1)
load(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_1/methData_mvalues_BIOS.RData")

#probe positions
row_gr <- rowRanges(mvalues_BIOS)

#split up by chromosome for easier processing
for (chr in 1:22) {
  chr_name <- paste0("chr", chr)

  #probes on this chromosome
  chr_idx <- seqnames(row_gr) == chr_name

  #subset the SummarizedExperiment
  mvalues_chr <- mvalues_BIOS[chr_idx, ]

  #report probe count
  message(chr_name, ": ", sum(chr_idx), " probes")

  #save per-chromosome object
  save_path <- paste0("/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_1/mvalues_BIOS_", chr_name, ".RData")
  save(mvalues_chr, file = save_path)
}

#load raw mRNA file 'counts'
load(file ="/exports/archive/kft-pgx-archive/BIOS/rnaSeqData_ReadCounts_BIOS_Freeze2_unrelated.RData")

#load PharmGKB gene list
genes <- read.delim("/home/dasmith/genes.tsv", stringsAsFactors = FALSE)

#parse to only annotated genes --> 2007 genes
genes <- genes[genes$Has.Variant.Annotation == "Yes", ]

#process mRNA counts: remove samples with missing flowcell
count.coldata <- colData(counts)
idx <- which(is.na(count.coldata$flowcell_num) == TRUE)
counts <- counts[,-idx]

#calculate RPKM values
counts_matrix <- assay(counts)
row_gr <- rowRanges(counts)
gene_lengths <- width(row_gr)
dge <- DGEList(counts = counts_matrix)
dge$genes$Length <- gene_lengths
rpkm_values <- rpkm(dge)
gene_avg_rpkm <- rowMeans(rpkm_values)

#filter to above 1 RPKM
filtered_gene_names <- names(gene_avg_rpkm[gene_avg_rpkm > 1])
counts_filtered <- counts[rownames(counts) %in% filtered_gene_names, ]
counts_pgkb <- counts_filtered[rownames(counts_filtered) %in% genes$Ensembl.Id, ] #only keep those annotated in PharmGKB (420 genes)

#replace raw counts with rpkm vals in filtered counts object
rpkm_pgkb <- rpkm_values[, colnames(counts_pgkb)]
rpkm_pgkb <- rpkm_pgkb[rownames(rpkm_pgkb) %in% filtered_gene_names, ]
rpkm_pgkb <- rpkm_pgkb[rownames(rpkm_pgkb) %in% genes$Ensembl.Id, ]
rpkm_se <- SummarizedExperiment(
  assays = list(rpkm = rpkm_pgkb),
  rowData = rowData(counts_pgkb),
  colData = colData(counts_pgkb),
  metadata = metadata(counts_pgkb)
)

#save mRNA file
save(rpkm_se, file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_4/rpkm_pgkb.RData")
