#title: DNA_processing
#author: Delaney Smith, adapted from Yunfeng
#date: '2025-2-26'
#function: Extract SNP dosages per cohort for chr1-22 (except chr16), merge cohorts, filter on MAF
#output: id lists (.txt), Output_02_mvalues_hrcv1.1.RData, chr<N>_BIOS_dosages.RData per chromosome

#load libraries
library(BiocParallel)
library(readr)
library(stringr)
library(dplyr)
library(SummarizedExperiment)

#set paths
path <- "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/"
path_geno<-"/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/" #add genotype file processing to this script

#load BIOS sample annotation (mvalues_BIOS) from script 1, 4031 people
load(file="/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_1/m_values_BIOS_chr16.RData")

#imputation ids (hrcv1.1) are not unique across the 6 cohorts, so extract per cohort
#per-cohort sample annotation: mvalues_hrcv1.1_col, mvalues_hrcv1.1
mvalues_col<-as.data.frame(colData(mvalues_BIOS))
mvalues_hrcv1.1_col<-mvalues_col[!is.na(mvalues_col$imputation_id_hrcv1.1),]
mvalues_hrcv1.1_col<-mvalues_hrcv1.1_col[,c("uuid","biobank_id","sampling_age","sex","gwas_chip","imputation_id_hrc","imputation_id_hrcv1.1")]
mvalues_hrcv1.1<-mvalues_BIOS[,rownames(mvalues_hrcv1.1_col)]

#fix LL cohort imputation id (raw id is wrong)
id_LL<-which(mvalues_hrcv1.1_col$biobank_id=="LL")
mvalues_hrcv1.1$imputation_id_hrcv1.1[id_LL]<-substring(mvalues_hrcv1.1$imputation_id_hrc[id_LL],10,13)

#remove duplicate and missing samples (gender mismatch), split by cohort and sex, write id lists
#produces mvalues_hrcv1.1_col_<cohort>[_male|_female]
mvalues_hrcv1.1_col_CODAM<-filter(mvalues_hrcv1.1_col,biobank_id=="CODAM")
mvalues_hrcv1.1_col_CODAM_male<-filter(mvalues_hrcv1.1_col_CODAM,sex=="male")
mvalues_hrcv1.1_col_CODAM_female<-filter(mvalues_hrcv1.1_col_CODAM,sex=="female")
write.table(mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_CODAM.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_CODAM_male$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_CODAM_male.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_CODAM_female$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_CODAM_female.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

mvalues_hrcv1.1_col_LL<-filter(mvalues_hrcv1.1_col,biobank_id=="LL")
mvalues_hrcv1.1_col_LL<-filter(mvalues_hrcv1.1_col_LL,!imputation_id_hrcv1.1%in%"1263") # 1263 sample is missing
mvalues_hrcv1.1_col_LL_male<-filter(mvalues_hrcv1.1_col_LL,sex=="male")
mvalues_hrcv1.1_col_LL_female<-filter(mvalues_hrcv1.1_col_LL,sex=="female")
write.table(mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_LL.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_LL_male$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_LL_male.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_LL_female$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_LL_female.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

mvalues_hrcv1.1_col_LLS<-filter(mvalues_hrcv1.1_col,biobank_id=="LLS")
mvalues_hrcv1.1_col_LLS<-filter(mvalues_hrcv1.1_col_LLS,!imputation_id_hrcv1.1%in%c("425","563","838","974","1436","2068")) # 6 samples are measured 2 times
mvalues_hrcv1.1_col_LLS_male<-filter(mvalues_hrcv1.1_col_LLS,sex=="male")
mvalues_hrcv1.1_col_LLS_female<-filter(mvalues_hrcv1.1_col_LLS,sex=="female")
write.table(mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_LLS.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_LLS_male$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_LLS_male.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_LLS_female$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_LLS_female.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

mvalues_hrcv1.1_col_NTR<-filter(mvalues_hrcv1.1_col,biobank_id=="NTR")
mvalues_hrcv1.1_col_NTR_male<-filter(mvalues_hrcv1.1_col_NTR,sex=="male")
mvalues_hrcv1.1_col_NTR_female<-filter(mvalues_hrcv1.1_col_NTR,sex=="female")
write.table(mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_NTR.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_NTR_male$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_NTR_male.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_NTR_female$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_NTR_female.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

mvalues_hrcv1.1_col_PAN<-filter(mvalues_hrcv1.1_col,biobank_id=="PAN")
mvalues_hrcv1.1_col_PAN_male<-filter(mvalues_hrcv1.1_col_PAN,sex=="male")
mvalues_hrcv1.1_col_PAN_female<-filter(mvalues_hrcv1.1_col_PAN,sex=="female")
write.table(mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_PAN.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_PAN_male$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_PAN_male.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_PAN_female$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_PAN_female.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

mvalues_hrcv1.1_col_RS<-filter(mvalues_hrcv1.1_col,biobank_id=="RS")
mvalues_hrcv1.1_col_RS_male<-filter(mvalues_hrcv1.1_col_RS,sex=="male")
mvalues_hrcv1.1_col_RS_female<-filter(mvalues_hrcv1.1_col_RS,sex=="female")
write.table(mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_RS.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_RS_male$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_RS_male.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mvalues_hrcv1.1_col_RS_female$imputation_id_hrcv1.1,file = "mvalues_hrcv1.1_ids_RS_female.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

#combine cohorts and subset mvalues_hrcv1.1
mvalues_hrcv1.1_col<-rbind(mvalues_hrcv1.1_col_CODAM,mvalues_hrcv1.1_col_LL,mvalues_hrcv1.1_col_LLS,mvalues_hrcv1.1_col_NTR,mvalues_hrcv1.1_col_PAN,mvalues_hrcv1.1_col_RS)
mvalues_hrcv1.1<-mvalues_hrcv1.1[,rownames(mvalues_hrcv1.1_col)]
save(mvalues_hrcv1.1,file = "Output_02_mvalues_hrcv1.1.RData") #already saved

#dosages per chromosome and cohort; objects dat_chr<N>_<cohort> and chr<N>_BIOS_dosages
#chr1: dat_chr1_<cohort> -> chr1_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr1_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr1.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr1_CODAM <- paste(dat_chr1_CODAM$CHROM,dat_chr1_CODAM$POS, sep = ":")
dat_chr1_CODAM$CHROM <- dat_chr1_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr1_CODAM),]
colnames(dat_chr1_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr1_CODAM <- as.matrix(dat_chr1_CODAM)
rownames(dat_chr1_CODAM) <- rowns_chr1_CODAM

#LL
dat_chr1_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr1.LL.DS.FORMAT", progress = FALSE)
rowns_chr1_LL <- paste(dat_chr1_LL$CHROM,dat_chr1_LL$POS, sep = ":")
dat_chr1_LL$CHROM <- dat_chr1_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr1_LL),]
colnames(dat_chr1_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr1_LL <- as.matrix(dat_chr1_LL)
rownames(dat_chr1_LL) <- rowns_chr1_LL

#PAN
dat_chr1_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr1.PAN.DS.FORMAT", progress = FALSE)
rowns_chr1_PAN <- paste(dat_chr1_PAN$CHROM, dat_chr1_PAN$POS, sep = ":")
dat_chr1_PAN$CHROM <- dat_chr1_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr1_PAN),]
colnames(dat_chr1_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr1_PAN <- as.matrix(dat_chr1_PAN)
rownames(dat_chr1_PAN) <- rowns_chr1_PAN

#LLS 660Q array
dat_chr1_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr1.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr1_LLS_660Q <- paste(dat_chr1_LLS_660Q$CHROM, dat_chr1_LLS_660Q$POS, sep = ":")
dat_chr1_LLS_660Q$CHROM <- dat_chr1_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr1_LLS_660Q),]
colnames(dat_chr1_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr1_LLS_660Q <- as.matrix(dat_chr1_LLS_660Q)
rownames(dat_chr1_LLS_660Q) <- rowns_chr1_LLS_660Q

#LLS OmniExpr array
dat_chr1_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr1.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr1_LLS_OmniExpr <- paste(dat_chr1_LLS_OmniExpr$CHROM, dat_chr1_LLS_OmniExpr$POS, sep = ":")
dat_chr1_LLS_OmniExpr$CHROM <- dat_chr1_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr1_LLS_OmniExpr),]
colnames(dat_chr1_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr1_LLS_OmniExpr <- as.matrix(dat_chr1_LLS_OmniExpr)
rownames(dat_chr1_LLS_OmniExpr) <- rowns_chr1_LLS_OmniExpr

#RS
dat_chr1_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr1.RS.DS.FORMAT", progress = FALSE)
rowns_chr1_RS <- paste(dat_chr1_RS$CHROM, dat_chr1_RS$POS, sep = ":")
dat_chr1_RS$CHROM <- dat_chr1_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr1_RS),]
colnames(dat_chr1_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr1_RS <- as.matrix(dat_chr1_RS)
rownames(dat_chr1_RS) <- rowns_chr1_RS

#NTR Affy6 array
dat_chr1_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr1.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr1_NTR_Affy6 <- paste(dat_chr1_NTR_Affy6$CHROM,dat_chr1_NTR_Affy6$POS, sep = ":")
dat_chr1_NTR_Affy6$CHROM <- dat_chr1_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr1_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr1_NTR_Affy6),]
colnames(dat_chr1_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr1_NTR_Affy6 <- as.matrix(dat_chr1_NTR_Affy6)
rownames(dat_chr1_NTR_Affy6) <- rowns_chr1_NTR_Affy6

#NTR GoNL array
dat_chr1_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr1.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr1_NTR_GoNL <- paste(dat_chr1_NTR_GoNL$CHROM, dat_chr1_NTR_GoNL$POS, sep = ":")
dat_chr1_NTR_GoNL$CHROM <- dat_chr1_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr1_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr1_NTR_GoNL),]
colnames(dat_chr1_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr1_NTR_GoNL <- as.matrix(dat_chr1_NTR_GoNL)
rownames(dat_chr1_NTR_GoNL) <- rowns_chr1_NTR_GoNL

#keep snps present in every cohort and bind samples
chr1_BIOS_dosages<-list(dat_chr1_CODAM,dat_chr1_LL,dat_chr1_LLS_660Q,dat_chr1_LLS_OmniExpr,dat_chr1_NTR_Affy6,dat_chr1_NTR_GoNL,dat_chr1_PAN,dat_chr1_RS)
rowns <- Reduce(intersect,lapply(chr1_BIOS_dosages, FUN = rownames))
chr1_BIOS_dosages <- lapply(chr1_BIOS_dosages, function(x) x[rowns,])
chr1_BIOS_dosages<-do.call(cbind,chr1_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr1_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr1_BIOS_dosages<-chr1_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr1_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr1_BIOS_dosages <- chr1_BIOS_dosages[, common_columns]
save(chr1_BIOS_dosages,file = "chr1_BIOS_dosages.RData") #yay done

#chr2: dat_chr2_<cohort> -> chr2_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr2_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr2.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr2_CODAM <- paste(dat_chr2_CODAM$CHROM,dat_chr2_CODAM$POS, sep = ":")
dat_chr2_CODAM$CHROM <- dat_chr2_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr2_CODAM),]
colnames(dat_chr2_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr2_CODAM <- as.matrix(dat_chr2_CODAM)
rownames(dat_chr2_CODAM) <- rowns_chr2_CODAM

#LL
dat_chr2_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr2.LL.DS.FORMAT", progress = FALSE)
rowns_chr2_LL <- paste(dat_chr2_LL$CHROM,dat_chr2_LL$POS, sep = ":")
dat_chr2_LL$CHROM <- dat_chr2_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr2_LL),]
colnames(dat_chr2_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr2_LL <- as.matrix(dat_chr2_LL)
rownames(dat_chr2_LL) <- rowns_chr2_LL

#PAN
dat_chr2_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr2.PAN.DS.FORMAT", progress = FALSE)
rowns_chr2_PAN <- paste(dat_chr2_PAN$CHROM, dat_chr2_PAN$POS, sep = ":")
dat_chr2_PAN$CHROM <- dat_chr2_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr2_PAN),]
colnames(dat_chr2_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr2_PAN <- as.matrix(dat_chr2_PAN)
rownames(dat_chr2_PAN) <- rowns_chr2_PAN

#LLS 660Q array
dat_chr2_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr2.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr2_LLS_660Q <- paste(dat_chr2_LLS_660Q$CHROM, dat_chr2_LLS_660Q$POS, sep = ":")
dat_chr2_LLS_660Q$CHROM <- dat_chr2_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr2_LLS_660Q),]
colnames(dat_chr2_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr2_LLS_660Q <- as.matrix(dat_chr2_LLS_660Q)
rownames(dat_chr2_LLS_660Q) <- rowns_chr2_LLS_660Q

#LLS OmniExpr array
dat_chr2_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr2.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr2_LLS_OmniExpr <- paste(dat_chr2_LLS_OmniExpr$CHROM, dat_chr2_LLS_OmniExpr$POS, sep = ":")
dat_chr2_LLS_OmniExpr$CHROM <- dat_chr2_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr2_LLS_OmniExpr),]
colnames(dat_chr2_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr2_LLS_OmniExpr <- as.matrix(dat_chr2_LLS_OmniExpr)
rownames(dat_chr2_LLS_OmniExpr) <- rowns_chr2_LLS_OmniExpr

#RS
dat_chr2_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr2.RS.DS.FORMAT", progress = FALSE)
rowns_chr2_RS <- paste(dat_chr2_RS$CHROM, dat_chr2_RS$POS, sep = ":")
dat_chr2_RS$CHROM <- dat_chr2_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr2_RS),]
colnames(dat_chr2_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr2_RS <- as.matrix(dat_chr2_RS)
rownames(dat_chr2_RS) <- rowns_chr2_RS

#NTR Affy6 array
dat_chr2_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr2.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr2_NTR_Affy6 <- paste(dat_chr2_NTR_Affy6$CHROM,dat_chr2_NTR_Affy6$POS, sep = ":")
dat_chr2_NTR_Affy6$CHROM <- dat_chr2_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr2_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr2_NTR_Affy6),]
colnames(dat_chr2_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr2_NTR_Affy6 <- as.matrix(dat_chr2_NTR_Affy6)
rownames(dat_chr2_NTR_Affy6) <- rowns_chr2_NTR_Affy6

#NTR GoNL array
dat_chr2_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr2.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr2_NTR_GoNL <- paste(dat_chr2_NTR_GoNL$CHROM, dat_chr2_NTR_GoNL$POS, sep = ":")
dat_chr2_NTR_GoNL$CHROM <- dat_chr2_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr2_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr2_NTR_GoNL),]
colnames(dat_chr2_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr2_NTR_GoNL <- as.matrix(dat_chr2_NTR_GoNL)
rownames(dat_chr2_NTR_GoNL) <- rowns_chr2_NTR_GoNL

#keep snps present in every cohort and bind samples
chr2_BIOS_dosages<-list(dat_chr2_CODAM,dat_chr2_LL,dat_chr2_LLS_660Q,dat_chr2_LLS_OmniExpr,dat_chr2_NTR_Affy6,dat_chr2_NTR_GoNL,dat_chr2_PAN,dat_chr2_RS)
rowns <- Reduce(intersect,lapply(chr2_BIOS_dosages, FUN = rownames))
chr2_BIOS_dosages <- lapply(chr2_BIOS_dosages, function(x) x[rowns,])
chr2_BIOS_dosages<-do.call(cbind,chr2_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr2_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr2_BIOS_dosages<-chr2_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr2_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr2_BIOS_dosages <- chr2_BIOS_dosages[, common_columns]
save(chr2_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr2_BIOS_dosages.RData") #yay done

#chr3: dat_chr3_<cohort> -> chr3_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr3_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr3.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr3_CODAM <- paste(dat_chr3_CODAM$CHROM,dat_chr3_CODAM$POS, sep = ":")
dat_chr3_CODAM$CHROM <- dat_chr3_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr3_CODAM),]
colnames(dat_chr3_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr3_CODAM <- as.matrix(dat_chr3_CODAM)
rownames(dat_chr3_CODAM) <- rowns_chr3_CODAM

#LL
dat_chr3_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr3.LL.DS.FORMAT", progress = FALSE)
rowns_chr3_LL <- paste(dat_chr3_LL$CHROM,dat_chr3_LL$POS, sep = ":")
dat_chr3_LL$CHROM <- dat_chr3_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr3_LL),]
colnames(dat_chr3_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr3_LL <- as.matrix(dat_chr3_LL)
rownames(dat_chr3_LL) <- rowns_chr3_LL

#PAN
dat_chr3_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr3.PAN.DS.FORMAT", progress = FALSE)
rowns_chr3_PAN <- paste(dat_chr3_PAN$CHROM, dat_chr3_PAN$POS, sep = ":")
dat_chr3_PAN$CHROM <- dat_chr3_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr3_PAN),]
colnames(dat_chr3_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr3_PAN <- as.matrix(dat_chr3_PAN)
rownames(dat_chr3_PAN) <- rowns_chr3_PAN

#LLS 660Q array
dat_chr3_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr3.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr3_LLS_660Q <- paste(dat_chr3_LLS_660Q$CHROM, dat_chr3_LLS_660Q$POS, sep = ":")
dat_chr3_LLS_660Q$CHROM <- dat_chr3_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr3_LLS_660Q),]
colnames(dat_chr3_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr3_LLS_660Q <- as.matrix(dat_chr3_LLS_660Q)
rownames(dat_chr3_LLS_660Q) <- rowns_chr3_LLS_660Q

#LLS OmniExpr array
dat_chr3_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr3.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr3_LLS_OmniExpr <- paste(dat_chr3_LLS_OmniExpr$CHROM, dat_chr3_LLS_OmniExpr$POS, sep = ":")
dat_chr3_LLS_OmniExpr$CHROM <- dat_chr3_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr3_LLS_OmniExpr),]
colnames(dat_chr3_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr3_LLS_OmniExpr <- as.matrix(dat_chr3_LLS_OmniExpr)
rownames(dat_chr3_LLS_OmniExpr) <- rowns_chr3_LLS_OmniExpr

#RS
dat_chr3_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr3.RS.DS.FORMAT", progress = FALSE)
rowns_chr3_RS <- paste(dat_chr3_RS$CHROM, dat_chr3_RS$POS, sep = ":")
dat_chr3_RS$CHROM <- dat_chr3_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr3_RS),]
colnames(dat_chr3_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr3_RS <- as.matrix(dat_chr3_RS)
rownames(dat_chr3_RS) <- rowns_chr3_RS

#NTR Affy6 array
dat_chr3_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr3.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr3_NTR_Affy6 <- paste(dat_chr3_NTR_Affy6$CHROM,dat_chr3_NTR_Affy6$POS, sep = ":")
dat_chr3_NTR_Affy6$CHROM <- dat_chr3_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr3_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr3_NTR_Affy6),]
colnames(dat_chr3_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr3_NTR_Affy6 <- as.matrix(dat_chr3_NTR_Affy6)
rownames(dat_chr3_NTR_Affy6) <- rowns_chr3_NTR_Affy6

#NTR GoNL array
dat_chr3_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr3.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr3_NTR_GoNL <- paste(dat_chr3_NTR_GoNL$CHROM, dat_chr3_NTR_GoNL$POS, sep = ":")
dat_chr3_NTR_GoNL$CHROM <- dat_chr3_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr3_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr3_NTR_GoNL),]
colnames(dat_chr3_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr3_NTR_GoNL <- as.matrix(dat_chr3_NTR_GoNL)
rownames(dat_chr3_NTR_GoNL) <- rowns_chr3_NTR_GoNL

#keep snps present in every cohort and bind samples
chr3_BIOS_dosages<-list(dat_chr3_CODAM,dat_chr3_LL,dat_chr3_LLS_660Q,dat_chr3_LLS_OmniExpr,dat_chr3_NTR_Affy6,dat_chr3_NTR_GoNL,dat_chr3_PAN,dat_chr3_RS)
rowns <- Reduce(intersect,lapply(chr3_BIOS_dosages, FUN = rownames))
chr3_BIOS_dosages <- lapply(chr3_BIOS_dosages, function(x) x[rowns,])
chr3_BIOS_dosages<-do.call(cbind,chr3_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr3_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr3_BIOS_dosages<-chr3_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr3_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr3_BIOS_dosages <- chr3_BIOS_dosages[, common_columns]
save(chr3_BIOS_dosages,file = "chr3_BIOS_dosages.RData") #yay done

#chr4: dat_chr4_<cohort> -> chr4_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr4_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr4.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr4_CODAM <- paste(dat_chr4_CODAM$CHROM,dat_chr4_CODAM$POS, sep = ":")
dat_chr4_CODAM$CHROM <- dat_chr4_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr4_CODAM),]
colnames(dat_chr4_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr4_CODAM <- as.matrix(dat_chr4_CODAM)
rownames(dat_chr4_CODAM) <- rowns_chr4_CODAM

#LL
dat_chr4_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr4.LL.DS.FORMAT", progress = FALSE)
rowns_chr4_LL <- paste(dat_chr4_LL$CHROM,dat_chr4_LL$POS, sep = ":")
dat_chr4_LL$CHROM <- dat_chr4_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr4_LL),]
colnames(dat_chr4_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr4_LL <- as.matrix(dat_chr4_LL)
rownames(dat_chr4_LL) <- rowns_chr4_LL

#PAN
dat_chr4_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr4.PAN.DS.FORMAT", progress = FALSE)
rowns_chr4_PAN <- paste(dat_chr4_PAN$CHROM, dat_chr4_PAN$POS, sep = ":")
dat_chr4_PAN$CHROM <- dat_chr4_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr4_PAN),]
colnames(dat_chr4_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr4_PAN <- as.matrix(dat_chr4_PAN)
rownames(dat_chr4_PAN) <- rowns_chr4_PAN

#LLS 660Q array
dat_chr4_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr4.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr4_LLS_660Q <- paste(dat_chr4_LLS_660Q$CHROM, dat_chr4_LLS_660Q$POS, sep = ":")
dat_chr4_LLS_660Q$CHROM <- dat_chr4_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr4_LLS_660Q),]
colnames(dat_chr4_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr4_LLS_660Q <- as.matrix(dat_chr4_LLS_660Q)
rownames(dat_chr4_LLS_660Q) <- rowns_chr4_LLS_660Q

#LLS OmniExpr array
dat_chr4_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr4.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr4_LLS_OmniExpr <- paste(dat_chr4_LLS_OmniExpr$CHROM, dat_chr4_LLS_OmniExpr$POS, sep = ":")
dat_chr4_LLS_OmniExpr$CHROM <- dat_chr4_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr4_LLS_OmniExpr),]
colnames(dat_chr4_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr4_LLS_OmniExpr <- as.matrix(dat_chr4_LLS_OmniExpr)
rownames(dat_chr4_LLS_OmniExpr) <- rowns_chr4_LLS_OmniExpr

#RS
dat_chr4_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr4.RS.DS.FORMAT", progress = FALSE)
rowns_chr4_RS <- paste(dat_chr4_RS$CHROM, dat_chr4_RS$POS, sep = ":")
dat_chr4_RS$CHROM <- dat_chr4_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr4_RS),]
colnames(dat_chr4_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr4_RS <- as.matrix(dat_chr4_RS)
rownames(dat_chr4_RS) <- rowns_chr4_RS

#NTR Affy6 array
dat_chr4_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr4.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr4_NTR_Affy6 <- paste(dat_chr4_NTR_Affy6$CHROM,dat_chr4_NTR_Affy6$POS, sep = ":")
dat_chr4_NTR_Affy6$CHROM <- dat_chr4_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr4_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr4_NTR_Affy6),]
colnames(dat_chr4_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr4_NTR_Affy6 <- as.matrix(dat_chr4_NTR_Affy6)
rownames(dat_chr4_NTR_Affy6) <- rowns_chr4_NTR_Affy6

#NTR GoNL array
dat_chr4_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr4.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr4_NTR_GoNL <- paste(dat_chr4_NTR_GoNL$CHROM, dat_chr4_NTR_GoNL$POS, sep = ":")
dat_chr4_NTR_GoNL$CHROM <- dat_chr4_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr4_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr4_NTR_GoNL),]
colnames(dat_chr4_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr4_NTR_GoNL <- as.matrix(dat_chr4_NTR_GoNL)
rownames(dat_chr4_NTR_GoNL) <- rowns_chr4_NTR_GoNL

#keep snps present in every cohort and bind samples
chr4_BIOS_dosages<-list(dat_chr4_CODAM,dat_chr4_LL,dat_chr4_LLS_660Q,dat_chr4_LLS_OmniExpr,dat_chr4_NTR_Affy6,dat_chr4_NTR_GoNL,dat_chr4_PAN,dat_chr4_RS)
rowns <- Reduce(intersect,lapply(chr4_BIOS_dosages, FUN = rownames))
chr4_BIOS_dosages <- lapply(chr4_BIOS_dosages, function(x) x[rowns,])
chr4_BIOS_dosages<-do.call(cbind,chr4_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr4_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr4_BIOS_dosages<-chr4_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr4_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr4_BIOS_dosages <- chr4_BIOS_dosages[, common_columns]
save(chr4_BIOS_dosages,file = "chr4_BIOS_dosages.RData") #yay done

#chr5: dat_chr5_<cohort> -> chr5_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr5_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr5.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr5_CODAM <- paste(dat_chr5_CODAM$CHROM,dat_chr5_CODAM$POS, sep = ":")
dat_chr5_CODAM$CHROM <- dat_chr5_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr5_CODAM),]
colnames(dat_chr5_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr5_CODAM <- as.matrix(dat_chr5_CODAM)
rownames(dat_chr5_CODAM) <- rowns_chr5_CODAM

#LL
dat_chr5_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr5.LL.DS.FORMAT", progress = FALSE)
rowns_chr5_LL <- paste(dat_chr5_LL$CHROM,dat_chr5_LL$POS, sep = ":")
dat_chr5_LL$CHROM <- dat_chr5_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr5_LL),]
colnames(dat_chr5_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr5_LL <- as.matrix(dat_chr5_LL)
rownames(dat_chr5_LL) <- rowns_chr5_LL

#PAN
dat_chr5_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr5.PAN.DS.FORMAT", progress = FALSE)
rowns_chr5_PAN <- paste(dat_chr5_PAN$CHROM, dat_chr5_PAN$POS, sep = ":")
dat_chr5_PAN$CHROM <- dat_chr5_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr5_PAN),]
colnames(dat_chr5_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr5_PAN <- as.matrix(dat_chr5_PAN)
rownames(dat_chr5_PAN) <- rowns_chr5_PAN

#LLS 660Q array
dat_chr5_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr5.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr5_LLS_660Q <- paste(dat_chr5_LLS_660Q$CHROM, dat_chr5_LLS_660Q$POS, sep = ":")
dat_chr5_LLS_660Q$CHROM <- dat_chr5_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr5_LLS_660Q),]
colnames(dat_chr5_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr5_LLS_660Q <- as.matrix(dat_chr5_LLS_660Q)
rownames(dat_chr5_LLS_660Q) <- rowns_chr5_LLS_660Q

#LLS OmniExpr array
dat_chr5_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr5.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr5_LLS_OmniExpr <- paste(dat_chr5_LLS_OmniExpr$CHROM, dat_chr5_LLS_OmniExpr$POS, sep = ":")
dat_chr5_LLS_OmniExpr$CHROM <- dat_chr5_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr5_LLS_OmniExpr),]
colnames(dat_chr5_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr5_LLS_OmniExpr <- as.matrix(dat_chr5_LLS_OmniExpr)
rownames(dat_chr5_LLS_OmniExpr) <- rowns_chr5_LLS_OmniExpr

#RS
dat_chr5_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr5.RS.DS.FORMAT", progress = FALSE)
rowns_chr5_RS <- paste(dat_chr5_RS$CHROM, dat_chr5_RS$POS, sep = ":")
dat_chr5_RS$CHROM <- dat_chr5_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr5_RS),]
colnames(dat_chr5_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr5_RS <- as.matrix(dat_chr5_RS)
rownames(dat_chr5_RS) <- rowns_chr5_RS

#NTR Affy6 array
dat_chr5_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr5.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr5_NTR_Affy6 <- paste(dat_chr5_NTR_Affy6$CHROM,dat_chr5_NTR_Affy6$POS, sep = ":")
dat_chr5_NTR_Affy6$CHROM <- dat_chr5_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr5_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr5_NTR_Affy6),]
colnames(dat_chr5_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr5_NTR_Affy6 <- as.matrix(dat_chr5_NTR_Affy6)
rownames(dat_chr5_NTR_Affy6) <- rowns_chr5_NTR_Affy6

#NTR GoNL array
dat_chr5_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr5.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr5_NTR_GoNL <- paste(dat_chr5_NTR_GoNL$CHROM, dat_chr5_NTR_GoNL$POS, sep = ":")
dat_chr5_NTR_GoNL$CHROM <- dat_chr5_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr5_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr5_NTR_GoNL),]
colnames(dat_chr5_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr5_NTR_GoNL <- as.matrix(dat_chr5_NTR_GoNL)
rownames(dat_chr5_NTR_GoNL) <- rowns_chr5_NTR_GoNL

#keep snps present in every cohort and bind samples
chr5_BIOS_dosages<-list(dat_chr5_CODAM,dat_chr5_LL,dat_chr5_LLS_660Q,dat_chr5_LLS_OmniExpr,dat_chr5_NTR_Affy6,dat_chr5_NTR_GoNL,dat_chr5_PAN,dat_chr5_RS)
rowns <- Reduce(intersect,lapply(chr5_BIOS_dosages, FUN = rownames))
chr5_BIOS_dosages <- lapply(chr5_BIOS_dosages, function(x) x[rowns,])
chr5_BIOS_dosages<-do.call(cbind,chr5_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr5_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr5_BIOS_dosages<-chr5_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr5_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr5_BIOS_dosages <- chr5_BIOS_dosages[, common_columns]
save(chr5_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr5_BIOS_dosages.RData") #yay done

#chr6: dat_chr6_<cohort> -> chr6_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr6_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr6.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr6_CODAM <- paste(dat_chr6_CODAM$CHROM,dat_chr6_CODAM$POS, sep = ":")
dat_chr6_CODAM$CHROM <- dat_chr6_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr6_CODAM),]
colnames(dat_chr6_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr6_CODAM <- as.matrix(dat_chr6_CODAM)
rownames(dat_chr6_CODAM) <- rowns_chr6_CODAM

#LL
dat_chr6_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr6.LL.DS.FORMAT", progress = FALSE)
rowns_chr6_LL <- paste(dat_chr6_LL$CHROM,dat_chr6_LL$POS, sep = ":")
dat_chr6_LL$CHROM <- dat_chr6_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr6_LL),]
colnames(dat_chr6_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr6_LL <- as.matrix(dat_chr6_LL)
rownames(dat_chr6_LL) <- rowns_chr6_LL

#PAN
dat_chr6_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr6.PAN.DS.FORMAT", progress = FALSE)
rowns_chr6_PAN <- paste(dat_chr6_PAN$CHROM, dat_chr6_PAN$POS, sep = ":")
dat_chr6_PAN$CHROM <- dat_chr6_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr6_PAN),]
colnames(dat_chr6_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr6_PAN <- as.matrix(dat_chr6_PAN)
rownames(dat_chr6_PAN) <- rowns_chr6_PAN

#LLS 660Q array
dat_chr6_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr6.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr6_LLS_660Q <- paste(dat_chr6_LLS_660Q$CHROM, dat_chr6_LLS_660Q$POS, sep = ":")
dat_chr6_LLS_660Q$CHROM <- dat_chr6_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr6_LLS_660Q),]
colnames(dat_chr6_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr6_LLS_660Q <- as.matrix(dat_chr6_LLS_660Q)
rownames(dat_chr6_LLS_660Q) <- rowns_chr6_LLS_660Q

#LLS OmniExpr array
dat_chr6_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr6.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr6_LLS_OmniExpr <- paste(dat_chr6_LLS_OmniExpr$CHROM, dat_chr6_LLS_OmniExpr$POS, sep = ":")
dat_chr6_LLS_OmniExpr$CHROM <- dat_chr6_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr6_LLS_OmniExpr),]
colnames(dat_chr6_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr6_LLS_OmniExpr <- as.matrix(dat_chr6_LLS_OmniExpr)
rownames(dat_chr6_LLS_OmniExpr) <- rowns_chr6_LLS_OmniExpr

#RS
dat_chr6_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr6.RS.DS.FORMAT", progress = FALSE)
rowns_chr6_RS <- paste(dat_chr6_RS$CHROM, dat_chr6_RS$POS, sep = ":")
dat_chr6_RS$CHROM <- dat_chr6_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr6_RS),]
colnames(dat_chr6_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr6_RS <- as.matrix(dat_chr6_RS)
rownames(dat_chr6_RS) <- rowns_chr6_RS

#NTR Affy6 array
dat_chr6_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr6.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr6_NTR_Affy6 <- paste(dat_chr6_NTR_Affy6$CHROM,dat_chr6_NTR_Affy6$POS, sep = ":")
dat_chr6_NTR_Affy6$CHROM <- dat_chr6_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr6_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr6_NTR_Affy6),]
colnames(dat_chr6_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr6_NTR_Affy6 <- as.matrix(dat_chr6_NTR_Affy6)
rownames(dat_chr6_NTR_Affy6) <- rowns_chr6_NTR_Affy6

#NTR GoNL array
dat_chr6_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr6.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr6_NTR_GoNL <- paste(dat_chr6_NTR_GoNL$CHROM, dat_chr6_NTR_GoNL$POS, sep = ":")
dat_chr6_NTR_GoNL$CHROM <- dat_chr6_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr6_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr6_NTR_GoNL),]
colnames(dat_chr6_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr6_NTR_GoNL <- as.matrix(dat_chr6_NTR_GoNL)
rownames(dat_chr6_NTR_GoNL) <- rowns_chr6_NTR_GoNL

#keep snps present in every cohort and bind samples
chr6_BIOS_dosages<-list(dat_chr6_CODAM,dat_chr6_LL,dat_chr6_LLS_660Q,dat_chr6_LLS_OmniExpr,dat_chr6_NTR_Affy6,dat_chr6_NTR_GoNL,dat_chr6_PAN,dat_chr6_RS)
rowns <- Reduce(intersect,lapply(chr6_BIOS_dosages, FUN = rownames))
chr6_BIOS_dosages <- lapply(chr6_BIOS_dosages, function(x) x[rowns,])
chr6_BIOS_dosages<-do.call(cbind,chr6_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr6_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr6_BIOS_dosages<-chr6_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr6_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr6_BIOS_dosages <- chr6_BIOS_dosages[, common_columns]
save(chr6_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr6_BIOS_dosages.RData") #yay done

#chr7: dat_chr7_<cohort> -> chr7_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr7_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr7.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr7_CODAM <- paste(dat_chr7_CODAM$CHROM,dat_chr7_CODAM$POS, sep = ":")
dat_chr7_CODAM$CHROM <- dat_chr7_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr7_CODAM),]
colnames(dat_chr7_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr7_CODAM <- as.matrix(dat_chr7_CODAM)
rownames(dat_chr7_CODAM) <- rowns_chr7_CODAM

#LL
dat_chr7_LL<-read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr7.LL.DS.FORMAT", progress = FALSE)
rowns_chr7_LL <- paste(dat_chr7_LL$CHROM,dat_chr7_LL$POS, sep = ":")
dat_chr7_LL$CHROM <- dat_chr7_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr7_LL),]
colnames(dat_chr7_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr7_LL <- as.matrix(dat_chr7_LL)
rownames(dat_chr7_LL) <- rowns_chr7_LL

#PAN
dat_chr7_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr7.PAN.DS.FORMAT", progress = FALSE)
rowns_chr7_PAN <- paste(dat_chr7_PAN$CHROM, dat_chr7_PAN$POS, sep = ":")
dat_chr7_PAN$CHROM <- dat_chr7_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr7_PAN),]
colnames(dat_chr7_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr7_PAN <- as.matrix(dat_chr7_PAN)
rownames(dat_chr7_PAN) <- rowns_chr7_PAN

#LLS 660Q array
dat_chr7_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr7.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr7_LLS_660Q <- paste(dat_chr7_LLS_660Q$CHROM, dat_chr7_LLS_660Q$POS, sep = ":")
dat_chr7_LLS_660Q$CHROM <- dat_chr7_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr7_LLS_660Q),]
colnames(dat_chr7_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr7_LLS_660Q <- as.matrix(dat_chr7_LLS_660Q)
rownames(dat_chr7_LLS_660Q) <- rowns_chr7_LLS_660Q

#LLS OmniExpr array
dat_chr7_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr7.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr7_LLS_OmniExpr <- paste(dat_chr7_LLS_OmniExpr$CHROM, dat_chr7_LLS_OmniExpr$POS, sep = ":")
dat_chr7_LLS_OmniExpr$CHROM <- dat_chr7_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr7_LLS_OmniExpr),]
colnames(dat_chr7_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr7_LLS_OmniExpr <- as.matrix(dat_chr7_LLS_OmniExpr)
rownames(dat_chr7_LLS_OmniExpr) <- rowns_chr7_LLS_OmniExpr

#RS
dat_chr7_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr7.RS.DS.FORMAT", progress = FALSE)
rowns_chr7_RS <- paste(dat_chr7_RS$CHROM, dat_chr7_RS$POS, sep = ":")
dat_chr7_RS$CHROM <- dat_chr7_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr7_RS),]
colnames(dat_chr7_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr7_RS <- as.matrix(dat_chr7_RS)
rownames(dat_chr7_RS) <- rowns_chr7_RS

#NTR Affy6 array
dat_chr7_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr7.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr7_NTR_Affy6 <- paste(dat_chr7_NTR_Affy6$CHROM,dat_chr7_NTR_Affy6$POS, sep = ":")
dat_chr7_NTR_Affy6$CHROM <- dat_chr7_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr7_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr7_NTR_Affy6),]
colnames(dat_chr7_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr7_NTR_Affy6 <- as.matrix(dat_chr7_NTR_Affy6)
rownames(dat_chr7_NTR_Affy6) <- rowns_chr7_NTR_Affy6

#NTR GoNL array
dat_chr7_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr7.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr7_NTR_GoNL <- paste(dat_chr7_NTR_GoNL$CHROM, dat_chr7_NTR_GoNL$POS, sep = ":")
dat_chr7_NTR_GoNL$CHROM <- dat_chr7_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr7_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr7_NTR_GoNL),]
colnames(dat_chr7_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr7_NTR_GoNL <- as.matrix(dat_chr7_NTR_GoNL)
rownames(dat_chr7_NTR_GoNL) <- rowns_chr7_NTR_GoNL

#keep snps present in every cohort and bind samples
chr7_BIOS_dosages<-list(dat_chr7_CODAM,dat_chr7_LL,dat_chr7_LLS_660Q,dat_chr7_LLS_OmniExpr,dat_chr7_NTR_Affy6,dat_chr7_NTR_GoNL,dat_chr7_PAN,dat_chr7_RS)
rowns <- Reduce(intersect,lapply(chr7_BIOS_dosages, FUN = rownames))
chr7_BIOS_dosages <- lapply(chr7_BIOS_dosages, function(x) x[rowns,])
chr7_BIOS_dosages<-do.call(cbind,chr7_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr7_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr7_BIOS_dosages<-chr7_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr7_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr7_BIOS_dosages <- chr7_BIOS_dosages[, common_columns]
save(chr7_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr7_BIOS_dosages.RData") #yay done

#chr8: dat_chr8_<cohort> -> chr8_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr8_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr8.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr8_CODAM <- paste(dat_chr8_CODAM$CHROM,dat_chr8_CODAM$POS, sep = ":")
dat_chr8_CODAM$CHROM <- dat_chr8_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr8_CODAM),]
colnames(dat_chr8_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr8_CODAM <- as.matrix(dat_chr8_CODAM)
rownames(dat_chr8_CODAM) <- rowns_chr8_CODAM

#LL
dat_chr8_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr8.LL.DS.FORMAT", progress = FALSE)
rowns_chr8_LL <- paste(dat_chr8_LL$CHROM,dat_chr8_LL$POS, sep = ":")
dat_chr8_LL$CHROM <- dat_chr8_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr8_LL),]
colnames(dat_chr8_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr8_LL <- as.matrix(dat_chr8_LL)
rownames(dat_chr8_LL) <- rowns_chr8_LL

#PAN
dat_chr8_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr8.PAN.DS.FORMAT", progress = FALSE)
rowns_chr8_PAN <- paste(dat_chr8_PAN$CHROM, dat_chr8_PAN$POS, sep = ":")
dat_chr8_PAN$CHROM <- dat_chr8_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr8_PAN),]
colnames(dat_chr8_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr8_PAN <- as.matrix(dat_chr8_PAN)
rownames(dat_chr8_PAN) <- rowns_chr8_PAN

#LLS 660Q array
dat_chr8_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr8.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr8_LLS_660Q <- paste(dat_chr8_LLS_660Q$CHROM, dat_chr8_LLS_660Q$POS, sep = ":")
dat_chr8_LLS_660Q$CHROM <- dat_chr8_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr8_LLS_660Q),]
colnames(dat_chr8_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr8_LLS_660Q <- as.matrix(dat_chr8_LLS_660Q)
rownames(dat_chr8_LLS_660Q) <- rowns_chr8_LLS_660Q

#LLS OmniExpr array
dat_chr8_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr8.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr8_LLS_OmniExpr <- paste(dat_chr8_LLS_OmniExpr$CHROM, dat_chr8_LLS_OmniExpr$POS, sep = ":")
dat_chr8_LLS_OmniExpr$CHROM <- dat_chr8_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr8_LLS_OmniExpr),]
colnames(dat_chr8_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr8_LLS_OmniExpr <- as.matrix(dat_chr8_LLS_OmniExpr)
rownames(dat_chr8_LLS_OmniExpr) <- rowns_chr8_LLS_OmniExpr

#RS
dat_chr8_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr8.RS.DS.FORMAT", progress = FALSE)
rowns_chr8_RS <- paste(dat_chr8_RS$CHROM, dat_chr8_RS$POS, sep = ":")
dat_chr8_RS$CHROM <- dat_chr8_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr8_RS),]
colnames(dat_chr8_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr8_RS <- as.matrix(dat_chr8_RS)
rownames(dat_chr8_RS) <- rowns_chr8_RS

#NTR Affy6 array
dat_chr8_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr8.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr8_NTR_Affy6 <- paste(dat_chr8_NTR_Affy6$CHROM,dat_chr8_NTR_Affy6$POS, sep = ":")
dat_chr8_NTR_Affy6$CHROM <- dat_chr8_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr8_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr8_NTR_Affy6),]
colnames(dat_chr8_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr8_NTR_Affy6 <- as.matrix(dat_chr8_NTR_Affy6)
rownames(dat_chr8_NTR_Affy6) <- rowns_chr8_NTR_Affy6

#NTR GoNL array
dat_chr8_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr8.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr8_NTR_GoNL <- paste(dat_chr8_NTR_GoNL$CHROM, dat_chr8_NTR_GoNL$POS, sep = ":")
dat_chr8_NTR_GoNL$CHROM <- dat_chr8_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr8_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr8_NTR_GoNL),]
colnames(dat_chr8_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr8_NTR_GoNL <- as.matrix(dat_chr8_NTR_GoNL)
rownames(dat_chr8_NTR_GoNL) <- rowns_chr8_NTR_GoNL

#keep snps present in every cohort and bind samples
chr8_BIOS_dosages<-list(dat_chr8_CODAM,dat_chr8_LL,dat_chr8_LLS_660Q,dat_chr8_LLS_OmniExpr,dat_chr8_NTR_Affy6,dat_chr8_NTR_GoNL,dat_chr8_PAN,dat_chr8_RS)
rowns <- Reduce(intersect,lapply(chr8_BIOS_dosages, FUN = rownames))
chr8_BIOS_dosages <- lapply(chr8_BIOS_dosages, function(x) x[rowns,])
chr8_BIOS_dosages<-do.call(cbind,chr8_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr8_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr8_BIOS_dosages<-chr8_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr8_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr8_BIOS_dosages <- chr8_BIOS_dosages[, common_columns]
save(chr8_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr8_BIOS_dosages.RData") #yay done

#chr9: dat_chr9_<cohort> -> chr9_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr9_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr9.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr9_CODAM <- paste(dat_chr9_CODAM$CHROM,dat_chr9_CODAM$POS, sep = ":")
dat_chr9_CODAM$CHROM <- dat_chr9_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr9_CODAM),]
colnames(dat_chr9_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr9_CODAM <- as.matrix(dat_chr9_CODAM)
rownames(dat_chr9_CODAM) <- rowns_chr9_CODAM

#LL
dat_chr9_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr9.LL.DS.FORMAT", progress = FALSE)
rowns_chr9_LL <- paste(dat_chr9_LL$CHROM,dat_chr9_LL$POS, sep = ":")
dat_chr9_LL$CHROM <- dat_chr9_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr9_LL),]
colnames(dat_chr9_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr9_LL <- as.matrix(dat_chr9_LL)
rownames(dat_chr9_LL) <- rowns_chr9_LL

#PAN
dat_chr9_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr9.PAN.DS.FORMAT", progress = FALSE)
rowns_chr9_PAN <- paste(dat_chr9_PAN$CHROM, dat_chr9_PAN$POS, sep = ":")
dat_chr9_PAN$CHROM <- dat_chr9_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr9_PAN),]
colnames(dat_chr9_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr9_PAN <- as.matrix(dat_chr9_PAN)
rownames(dat_chr9_PAN) <- rowns_chr9_PAN

#LLS 660Q array
dat_chr9_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr9.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr9_LLS_660Q <- paste(dat_chr9_LLS_660Q$CHROM, dat_chr9_LLS_660Q$POS, sep = ":")
dat_chr9_LLS_660Q$CHROM <- dat_chr9_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr9_LLS_660Q),]
colnames(dat_chr9_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr9_LLS_660Q <- as.matrix(dat_chr9_LLS_660Q)
rownames(dat_chr9_LLS_660Q) <- rowns_chr9_LLS_660Q

#LLS OmniExpr array
dat_chr9_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr9.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr9_LLS_OmniExpr <- paste(dat_chr9_LLS_OmniExpr$CHROM, dat_chr9_LLS_OmniExpr$POS, sep = ":")
dat_chr9_LLS_OmniExpr$CHROM <- dat_chr9_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr9_LLS_OmniExpr),]
colnames(dat_chr9_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr9_LLS_OmniExpr <- as.matrix(dat_chr9_LLS_OmniExpr)
rownames(dat_chr9_LLS_OmniExpr) <- rowns_chr9_LLS_OmniExpr

#RS
dat_chr9_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr9.RS.DS.FORMAT", progress = FALSE)
rowns_chr9_RS <- paste(dat_chr9_RS$CHROM, dat_chr9_RS$POS, sep = ":")
dat_chr9_RS$CHROM <- dat_chr9_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr9_RS),]
colnames(dat_chr9_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr9_RS <- as.matrix(dat_chr9_RS)
rownames(dat_chr9_RS) <- rowns_chr9_RS

#NTR Affy6 array
dat_chr9_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr9.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr9_NTR_Affy6 <- paste(dat_chr9_NTR_Affy6$CHROM,dat_chr9_NTR_Affy6$POS, sep = ":")
dat_chr9_NTR_Affy6$CHROM <- dat_chr9_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr9_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr9_NTR_Affy6),]
colnames(dat_chr9_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr9_NTR_Affy6 <- as.matrix(dat_chr9_NTR_Affy6)
rownames(dat_chr9_NTR_Affy6) <- rowns_chr9_NTR_Affy6

#NTR GoNL array
dat_chr9_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr9.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr9_NTR_GoNL <- paste(dat_chr9_NTR_GoNL$CHROM, dat_chr9_NTR_GoNL$POS, sep = ":")
dat_chr9_NTR_GoNL$CHROM <- dat_chr9_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr9_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr9_NTR_GoNL),]
colnames(dat_chr9_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr9_NTR_GoNL <- as.matrix(dat_chr9_NTR_GoNL)
rownames(dat_chr9_NTR_GoNL) <- rowns_chr9_NTR_GoNL

#keep snps present in every cohort and bind samples
chr9_BIOS_dosages<-list(dat_chr9_CODAM,dat_chr9_LL,dat_chr9_LLS_660Q,dat_chr9_LLS_OmniExpr,dat_chr9_NTR_Affy6,dat_chr9_NTR_GoNL,dat_chr9_PAN,dat_chr9_RS)
rowns <- Reduce(intersect,lapply(chr9_BIOS_dosages, FUN = rownames))
chr9_BIOS_dosages <- lapply(chr9_BIOS_dosages, function(x) x[rowns,])
chr9_BIOS_dosages<-do.call(cbind,chr9_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr9_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr9_BIOS_dosages<-chr9_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr9_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr9_BIOS_dosages <- chr9_BIOS_dosages[, common_columns]
save(chr9_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr9_BIOS_dosages.RData") #yay done

#chr10: dat_chr10_<cohort> -> chr10_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr10_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr10.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr10_CODAM <- paste(dat_chr10_CODAM$CHROM,dat_chr10_CODAM$POS, sep = ":")
dat_chr10_CODAM$CHROM <- dat_chr10_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr10_CODAM),]
colnames(dat_chr10_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr10_CODAM <- as.matrix(dat_chr10_CODAM)
rownames(dat_chr10_CODAM) <- rowns_chr10_CODAM

#LL
dat_chr10_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr10.LL.DS.FORMAT", progress = FALSE)
rowns_chr10_LL <- paste(dat_chr10_LL$CHROM,dat_chr10_LL$POS, sep = ":")
dat_chr10_LL$CHROM <- dat_chr10_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr10_LL),]
colnames(dat_chr10_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr10_LL <- as.matrix(dat_chr10_LL)
rownames(dat_chr10_LL) <- rowns_chr10_LL

#PAN
dat_chr10_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr10.PAN.DS.FORMAT", progress = FALSE)
rowns_chr10_PAN <- paste(dat_chr10_PAN$CHROM, dat_chr10_PAN$POS, sep = ":")
dat_chr10_PAN$CHROM <- dat_chr10_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr10_PAN),]
colnames(dat_chr10_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr10_PAN <- as.matrix(dat_chr10_PAN)
rownames(dat_chr10_PAN) <- rowns_chr10_PAN

#LLS 660Q array
dat_chr10_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr10.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr10_LLS_660Q <- paste(dat_chr10_LLS_660Q$CHROM, dat_chr10_LLS_660Q$POS, sep = ":")
dat_chr10_LLS_660Q$CHROM <- dat_chr10_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr10_LLS_660Q),]
colnames(dat_chr10_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr10_LLS_660Q <- as.matrix(dat_chr10_LLS_660Q)
rownames(dat_chr10_LLS_660Q) <- rowns_chr10_LLS_660Q

#LLS OmniExpr array
dat_chr10_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr10.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr10_LLS_OmniExpr <- paste(dat_chr10_LLS_OmniExpr$CHROM, dat_chr10_LLS_OmniExpr$POS, sep = ":")
dat_chr10_LLS_OmniExpr$CHROM <- dat_chr10_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr10_LLS_OmniExpr),]
colnames(dat_chr10_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr10_LLS_OmniExpr <- as.matrix(dat_chr10_LLS_OmniExpr)
rownames(dat_chr10_LLS_OmniExpr) <- rowns_chr10_LLS_OmniExpr

#RS
dat_chr10_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr10.RS.DS.FORMAT", progress = FALSE)
rowns_chr10_RS <- paste(dat_chr10_RS$CHROM, dat_chr10_RS$POS, sep = ":")
dat_chr10_RS$CHROM <- dat_chr10_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr10_RS),]
colnames(dat_chr10_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr10_RS <- as.matrix(dat_chr10_RS)
rownames(dat_chr10_RS) <- rowns_chr10_RS

#NTR Affy6 array
dat_chr10_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr10.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr10_NTR_Affy6 <- paste(dat_chr10_NTR_Affy6$CHROM,dat_chr10_NTR_Affy6$POS, sep = ":")
dat_chr10_NTR_Affy6$CHROM <- dat_chr10_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr10_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr10_NTR_Affy6),]
colnames(dat_chr10_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr10_NTR_Affy6 <- as.matrix(dat_chr10_NTR_Affy6)
rownames(dat_chr10_NTR_Affy6) <- rowns_chr10_NTR_Affy6

#NTR GoNL array
dat_chr10_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr10.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr10_NTR_GoNL <- paste(dat_chr10_NTR_GoNL$CHROM, dat_chr10_NTR_GoNL$POS, sep = ":")
dat_chr10_NTR_GoNL$CHROM <- dat_chr10_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr10_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr10_NTR_GoNL),]
colnames(dat_chr10_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr10_NTR_GoNL <- as.matrix(dat_chr10_NTR_GoNL)
rownames(dat_chr10_NTR_GoNL) <- rowns_chr10_NTR_GoNL

#keep snps present in every cohort and bind samples
chr10_BIOS_dosages<-list(dat_chr10_CODAM,dat_chr10_LL,dat_chr10_LLS_660Q,dat_chr10_LLS_OmniExpr,dat_chr10_NTR_Affy6,dat_chr10_NTR_GoNL,dat_chr10_PAN,dat_chr10_RS)
rowns <- Reduce(intersect,lapply(chr10_BIOS_dosages, FUN = rownames))
chr10_BIOS_dosages <- lapply(chr10_BIOS_dosages, function(x) x[rowns,])
chr10_BIOS_dosages<-do.call(cbind,chr10_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr10_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr10_BIOS_dosages<-chr10_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr10_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr10_BIOS_dosages <- chr10_BIOS_dosages[, common_columns]
save(chr10_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr10_BIOS_dosages.RData") #yay done

#chr11: dat_chr11_<cohort> -> chr11_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr11_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr11.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr11_CODAM <- paste(dat_chr11_CODAM$CHROM,dat_chr11_CODAM$POS, sep = ":")
dat_chr11_CODAM$CHROM <- dat_chr11_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr11_CODAM),]
colnames(dat_chr11_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr11_CODAM <- as.matrix(dat_chr11_CODAM)
rownames(dat_chr11_CODAM) <- rowns_chr11_CODAM

#LL
dat_chr11_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr11.LL.DS.FORMAT", progress = FALSE)
rowns_chr11_LL <- paste(dat_chr11_LL$CHROM,dat_chr11_LL$POS, sep = ":")
dat_chr11_LL$CHROM <- dat_chr11_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr11_LL),]
colnames(dat_chr11_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr11_LL <- as.matrix(dat_chr11_LL)
rownames(dat_chr11_LL) <- rowns_chr11_LL

#PAN
dat_chr11_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr11.PAN.DS.FORMAT", progress = FALSE)
rowns_chr11_PAN <- paste(dat_chr11_PAN$CHROM, dat_chr11_PAN$POS, sep = ":")
dat_chr11_PAN$CHROM <- dat_chr11_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr11_PAN),]
colnames(dat_chr11_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr11_PAN <- as.matrix(dat_chr11_PAN)
rownames(dat_chr11_PAN) <- rowns_chr11_PAN

#LLS 660Q array
dat_chr11_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr11.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr11_LLS_660Q <- paste(dat_chr11_LLS_660Q$CHROM, dat_chr11_LLS_660Q$POS, sep = ":")
dat_chr11_LLS_660Q$CHROM <- dat_chr11_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr11_LLS_660Q),]
colnames(dat_chr11_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr11_LLS_660Q <- as.matrix(dat_chr11_LLS_660Q)
rownames(dat_chr11_LLS_660Q) <- rowns_chr11_LLS_660Q

#LLS OmniExpr array
dat_chr11_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr11.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr11_LLS_OmniExpr <- paste(dat_chr11_LLS_OmniExpr$CHROM, dat_chr11_LLS_OmniExpr$POS, sep = ":")
dat_chr11_LLS_OmniExpr$CHROM <- dat_chr11_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr11_LLS_OmniExpr),]
colnames(dat_chr11_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr11_LLS_OmniExpr <- as.matrix(dat_chr11_LLS_OmniExpr)
rownames(dat_chr11_LLS_OmniExpr) <- rowns_chr11_LLS_OmniExpr

#RS
dat_chr11_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr11.RS.DS.FORMAT", progress = FALSE)
rowns_chr11_RS <- paste(dat_chr11_RS$CHROM, dat_chr11_RS$POS, sep = ":")
dat_chr11_RS$CHROM <- dat_chr11_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr11_RS),]
colnames(dat_chr11_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr11_RS <- as.matrix(dat_chr11_RS)
rownames(dat_chr11_RS) <- rowns_chr11_RS

#NTR Affy6 array
dat_chr11_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr11.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr11_NTR_Affy6 <- paste(dat_chr11_NTR_Affy6$CHROM,dat_chr11_NTR_Affy6$POS, sep = ":")
dat_chr11_NTR_Affy6$CHROM <- dat_chr11_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr11_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr11_NTR_Affy6),]
colnames(dat_chr11_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr11_NTR_Affy6 <- as.matrix(dat_chr11_NTR_Affy6)
rownames(dat_chr11_NTR_Affy6) <- rowns_chr11_NTR_Affy6

#NTR GoNL array
dat_chr11_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr11.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr11_NTR_GoNL <- paste(dat_chr11_NTR_GoNL$CHROM, dat_chr11_NTR_GoNL$POS, sep = ":")
dat_chr11_NTR_GoNL$CHROM <- dat_chr11_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr11_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr11_NTR_GoNL),]
colnames(dat_chr11_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr11_NTR_GoNL <- as.matrix(dat_chr11_NTR_GoNL)
rownames(dat_chr11_NTR_GoNL) <- rowns_chr11_NTR_GoNL

#keep snps present in every cohort and bind samples
chr11_BIOS_dosages<-list(dat_chr11_CODAM,dat_chr11_LL,dat_chr11_LLS_660Q,dat_chr11_LLS_OmniExpr,dat_chr11_NTR_Affy6,dat_chr11_NTR_GoNL,dat_chr11_PAN,dat_chr11_RS)
rowns <- Reduce(intersect,lapply(chr11_BIOS_dosages, FUN = rownames))
chr11_BIOS_dosages <- lapply(chr11_BIOS_dosages, function(x) x[rowns,])
chr11_BIOS_dosages<-do.call(cbind,chr11_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr11_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr11_BIOS_dosages<-chr11_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr11_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr11_BIOS_dosages <- chr11_BIOS_dosages[, common_columns]
save(chr11_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr11_BIOS_dosages.RData") #yay done

#chr12: dat_chr12_<cohort> -> chr12_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr12_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr12.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr12_CODAM <- paste(dat_chr12_CODAM$CHROM,dat_chr12_CODAM$POS, sep = ":")
dat_chr12_CODAM$CHROM <- dat_chr12_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr12_CODAM),]
colnames(dat_chr12_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr12_CODAM <- as.matrix(dat_chr12_CODAM)
rownames(dat_chr12_CODAM) <- rowns_chr12_CODAM

#LL
dat_chr12_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr12.LL.DS.FORMAT", progress = FALSE)
rowns_chr12_LL <- paste(dat_chr12_LL$CHROM,dat_chr12_LL$POS, sep = ":")
dat_chr12_LL$CHROM <- dat_chr12_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr12_LL),]
colnames(dat_chr12_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr12_LL <- as.matrix(dat_chr12_LL)
rownames(dat_chr12_LL) <- rowns_chr12_LL

#PAN
dat_chr12_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr12.PAN.DS.FORMAT", progress = FALSE)
rowns_chr12_PAN <- paste(dat_chr12_PAN$CHROM, dat_chr12_PAN$POS, sep = ":")
dat_chr12_PAN$CHROM <- dat_chr12_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr12_PAN),]
colnames(dat_chr12_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr12_PAN <- as.matrix(dat_chr12_PAN)
rownames(dat_chr12_PAN) <- rowns_chr12_PAN

#LLS 660Q array
dat_chr12_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr12.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr12_LLS_660Q <- paste(dat_chr12_LLS_660Q$CHROM, dat_chr12_LLS_660Q$POS, sep = ":")
dat_chr12_LLS_660Q$CHROM <- dat_chr12_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr12_LLS_660Q),]
colnames(dat_chr12_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr12_LLS_660Q <- as.matrix(dat_chr12_LLS_660Q)
rownames(dat_chr12_LLS_660Q) <- rowns_chr12_LLS_660Q

#LLS OmniExpr array
dat_chr12_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr12.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr12_LLS_OmniExpr <- paste(dat_chr12_LLS_OmniExpr$CHROM, dat_chr12_LLS_OmniExpr$POS, sep = ":")
dat_chr12_LLS_OmniExpr$CHROM <- dat_chr12_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr12_LLS_OmniExpr),]
colnames(dat_chr12_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr12_LLS_OmniExpr <- as.matrix(dat_chr12_LLS_OmniExpr)
rownames(dat_chr12_LLS_OmniExpr) <- rowns_chr12_LLS_OmniExpr

#RS
dat_chr12_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr12.RS.DS.FORMAT", progress = FALSE)
rowns_chr12_RS <- paste(dat_chr12_RS$CHROM, dat_chr12_RS$POS, sep = ":")
dat_chr12_RS$CHROM <- dat_chr12_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr12_RS),]
colnames(dat_chr12_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr12_RS <- as.matrix(dat_chr12_RS)
rownames(dat_chr12_RS) <- rowns_chr12_RS

#NTR Affy6 array
dat_chr12_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr12.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr12_NTR_Affy6 <- paste(dat_chr12_NTR_Affy6$CHROM,dat_chr12_NTR_Affy6$POS, sep = ":")
dat_chr12_NTR_Affy6$CHROM <- dat_chr12_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr12_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr12_NTR_Affy6),]
colnames(dat_chr12_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr12_NTR_Affy6 <- as.matrix(dat_chr12_NTR_Affy6)
rownames(dat_chr12_NTR_Affy6) <- rowns_chr12_NTR_Affy6

#NTR GoNL array
dat_chr12_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr12.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr12_NTR_GoNL <- paste(dat_chr12_NTR_GoNL$CHROM, dat_chr12_NTR_GoNL$POS, sep = ":")
dat_chr12_NTR_GoNL$CHROM <- dat_chr12_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr12_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr12_NTR_GoNL),]
colnames(dat_chr12_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr12_NTR_GoNL <- as.matrix(dat_chr12_NTR_GoNL)
rownames(dat_chr12_NTR_GoNL) <- rowns_chr12_NTR_GoNL

#keep snps present in every cohort and bind samples
chr12_BIOS_dosages<-list(dat_chr12_CODAM,dat_chr12_LL,dat_chr12_LLS_660Q,dat_chr12_LLS_OmniExpr,dat_chr12_NTR_Affy6,dat_chr12_NTR_GoNL,dat_chr12_PAN,dat_chr12_RS)
rowns <- Reduce(intersect,lapply(chr12_BIOS_dosages, FUN = rownames))
chr12_BIOS_dosages <- lapply(chr12_BIOS_dosages, function(x) x[rowns,])
chr12_BIOS_dosages<-do.call(cbind,chr12_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr12_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr12_BIOS_dosages<-chr12_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr12_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr12_BIOS_dosages <- chr12_BIOS_dosages[, common_columns]
save(chr12_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr12_BIOS_dosages.RData") #yay done

#chr13: dat_chr13_<cohort> -> chr13_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr13_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr13.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr13_CODAM <- paste(dat_chr13_CODAM$CHROM,dat_chr13_CODAM$POS, sep = ":")
dat_chr13_CODAM$CHROM <- dat_chr13_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr13_CODAM),]
colnames(dat_chr13_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr13_CODAM <- as.matrix(dat_chr13_CODAM)
rownames(dat_chr13_CODAM) <- rowns_chr13_CODAM

#LL
dat_chr13_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr13.LL.DS.FORMAT", progress = FALSE)
rowns_chr13_LL <- paste(dat_chr13_LL$CHROM,dat_chr13_LL$POS, sep = ":")
dat_chr13_LL$CHROM <- dat_chr13_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr13_LL),]
colnames(dat_chr13_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr13_LL <- as.matrix(dat_chr13_LL)
rownames(dat_chr13_LL) <- rowns_chr13_LL

#PAN
dat_chr13_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr13.PAN.DS.FORMAT", progress = FALSE)
rowns_chr13_PAN <- paste(dat_chr13_PAN$CHROM, dat_chr13_PAN$POS, sep = ":")
dat_chr13_PAN$CHROM <- dat_chr13_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr13_PAN),]
colnames(dat_chr13_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr13_PAN <- as.matrix(dat_chr13_PAN)
rownames(dat_chr13_PAN) <- rowns_chr13_PAN

#LLS 660Q array
dat_chr13_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr13.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr13_LLS_660Q <- paste(dat_chr13_LLS_660Q$CHROM, dat_chr13_LLS_660Q$POS, sep = ":")
dat_chr13_LLS_660Q$CHROM <- dat_chr13_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr13_LLS_660Q),]
colnames(dat_chr13_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr13_LLS_660Q <- as.matrix(dat_chr13_LLS_660Q)
rownames(dat_chr13_LLS_660Q) <- rowns_chr13_LLS_660Q

#LLS OmniExpr array
dat_chr13_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr13.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr13_LLS_OmniExpr <- paste(dat_chr13_LLS_OmniExpr$CHROM, dat_chr13_LLS_OmniExpr$POS, sep = ":")
dat_chr13_LLS_OmniExpr$CHROM <- dat_chr13_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr13_LLS_OmniExpr),]
colnames(dat_chr13_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr13_LLS_OmniExpr <- as.matrix(dat_chr13_LLS_OmniExpr)
rownames(dat_chr13_LLS_OmniExpr) <- rowns_chr13_LLS_OmniExpr

#RS
dat_chr13_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr13.RS.DS.FORMAT", progress = FALSE)
rowns_chr13_RS <- paste(dat_chr13_RS$CHROM, dat_chr13_RS$POS, sep = ":")
dat_chr13_RS$CHROM <- dat_chr13_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr13_RS),]
colnames(dat_chr13_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr13_RS <- as.matrix(dat_chr13_RS)
rownames(dat_chr13_RS) <- rowns_chr13_RS

#NTR Affy6 array
dat_chr13_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr13.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr13_NTR_Affy6 <- paste(dat_chr13_NTR_Affy6$CHROM,dat_chr13_NTR_Affy6$POS, sep = ":")
dat_chr13_NTR_Affy6$CHROM <- dat_chr13_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr13_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr13_NTR_Affy6),]
colnames(dat_chr13_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr13_NTR_Affy6 <- as.matrix(dat_chr13_NTR_Affy6)
rownames(dat_chr13_NTR_Affy6) <- rowns_chr13_NTR_Affy6

#NTR GoNL array
dat_chr13_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr13.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr13_NTR_GoNL <- paste(dat_chr13_NTR_GoNL$CHROM, dat_chr13_NTR_GoNL$POS, sep = ":")
dat_chr13_NTR_GoNL$CHROM <- dat_chr13_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr13_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr13_NTR_GoNL),]
colnames(dat_chr13_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr13_NTR_GoNL <- as.matrix(dat_chr13_NTR_GoNL)
rownames(dat_chr13_NTR_GoNL) <- rowns_chr13_NTR_GoNL

#keep snps present in every cohort and bind samples
chr13_BIOS_dosages<-list(dat_chr13_CODAM,dat_chr13_LL,dat_chr13_LLS_660Q,dat_chr13_LLS_OmniExpr,dat_chr13_NTR_Affy6,dat_chr13_NTR_GoNL,dat_chr13_PAN,dat_chr13_RS)
rowns <- Reduce(intersect,lapply(chr13_BIOS_dosages, FUN = rownames))
chr13_BIOS_dosages <- lapply(chr13_BIOS_dosages, function(x) x[rowns,])
chr13_BIOS_dosages<-do.call(cbind,chr13_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr13_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr13_BIOS_dosages<-chr13_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr13_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr13_BIOS_dosages <- chr13_BIOS_dosages[, common_columns]
save(chr13_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr13_BIOS_dosages.RData") #yay done

#chr14: dat_chr14_<cohort> -> chr14_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr14_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr14.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr14_CODAM <- paste(dat_chr14_CODAM$CHROM,dat_chr14_CODAM$POS, sep = ":")
dat_chr14_CODAM$CHROM <- dat_chr14_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr14_CODAM),]
colnames(dat_chr14_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr14_CODAM <- as.matrix(dat_chr14_CODAM)
rownames(dat_chr14_CODAM) <- rowns_chr14_CODAM

#LL
dat_chr14_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr14.LL.DS.FORMAT", progress = FALSE)
rowns_chr14_LL <- paste(dat_chr14_LL$CHROM,dat_chr14_LL$POS, sep = ":")
dat_chr14_LL$CHROM <- dat_chr14_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr14_LL),]
colnames(dat_chr14_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr14_LL <- as.matrix(dat_chr14_LL)
rownames(dat_chr14_LL) <- rowns_chr14_LL

#PAN
dat_chr14_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr14.PAN.DS.FORMAT", progress = FALSE)
rowns_chr14_PAN <- paste(dat_chr14_PAN$CHROM, dat_chr14_PAN$POS, sep = ":")
dat_chr14_PAN$CHROM <- dat_chr14_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr14_PAN),]
colnames(dat_chr14_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr14_PAN <- as.matrix(dat_chr14_PAN)
rownames(dat_chr14_PAN) <- rowns_chr14_PAN

#LLS 660Q array
dat_chr14_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr14.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr14_LLS_660Q <- paste(dat_chr14_LLS_660Q$CHROM, dat_chr14_LLS_660Q$POS, sep = ":")
dat_chr14_LLS_660Q$CHROM <- dat_chr14_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr14_LLS_660Q),]
colnames(dat_chr14_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr14_LLS_660Q <- as.matrix(dat_chr14_LLS_660Q)
rownames(dat_chr14_LLS_660Q) <- rowns_chr14_LLS_660Q

#LLS OmniExpr array
dat_chr14_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr14.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr14_LLS_OmniExpr <- paste(dat_chr14_LLS_OmniExpr$CHROM, dat_chr14_LLS_OmniExpr$POS, sep = ":")
dat_chr14_LLS_OmniExpr$CHROM <- dat_chr14_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr14_LLS_OmniExpr),]
colnames(dat_chr14_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr14_LLS_OmniExpr <- as.matrix(dat_chr14_LLS_OmniExpr)
rownames(dat_chr14_LLS_OmniExpr) <- rowns_chr14_LLS_OmniExpr

#RS
dat_chr14_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr14.RS.DS.FORMAT", progress = FALSE)
rowns_chr14_RS <- paste(dat_chr14_RS$CHROM, dat_chr14_RS$POS, sep = ":")
dat_chr14_RS$CHROM <- dat_chr14_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr14_RS),]
colnames(dat_chr14_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr14_RS <- as.matrix(dat_chr14_RS)
rownames(dat_chr14_RS) <- rowns_chr14_RS

#NTR Affy6 array
dat_chr14_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr14.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr14_NTR_Affy6 <- paste(dat_chr14_NTR_Affy6$CHROM,dat_chr14_NTR_Affy6$POS, sep = ":")
dat_chr14_NTR_Affy6$CHROM <- dat_chr14_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr14_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr14_NTR_Affy6),]
colnames(dat_chr14_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr14_NTR_Affy6 <- as.matrix(dat_chr14_NTR_Affy6)
rownames(dat_chr14_NTR_Affy6) <- rowns_chr14_NTR_Affy6

#NTR GoNL array
dat_chr14_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr14.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr14_NTR_GoNL <- paste(dat_chr14_NTR_GoNL$CHROM, dat_chr14_NTR_GoNL$POS, sep = ":")
dat_chr14_NTR_GoNL$CHROM <- dat_chr14_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr14_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr14_NTR_GoNL),]
colnames(dat_chr14_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr14_NTR_GoNL <- as.matrix(dat_chr14_NTR_GoNL)
rownames(dat_chr14_NTR_GoNL) <- rowns_chr14_NTR_GoNL

#keep snps present in every cohort and bind samples
chr14_BIOS_dosages<-list(dat_chr14_CODAM,dat_chr14_LL,dat_chr14_LLS_660Q,dat_chr14_LLS_OmniExpr,dat_chr14_NTR_Affy6,dat_chr14_NTR_GoNL,dat_chr14_PAN,dat_chr14_RS)
rowns <- Reduce(intersect,lapply(chr14_BIOS_dosages, FUN = rownames))
chr14_BIOS_dosages <- lapply(chr14_BIOS_dosages, function(x) x[rowns,])
chr14_BIOS_dosages<-do.call(cbind,chr14_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr14_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr14_BIOS_dosages<-chr14_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr14_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr14_BIOS_dosages <- chr14_BIOS_dosages[, common_columns]
save(chr14_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr14_BIOS_dosages.RData") #yay done

#chr15: dat_chr15_<cohort> -> chr15_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr15_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr15.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr15_CODAM <- paste(dat_chr15_CODAM$CHROM,dat_chr15_CODAM$POS, sep = ":")
dat_chr15_CODAM$CHROM <- dat_chr15_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr15_CODAM),]
colnames(dat_chr15_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr15_CODAM <- as.matrix(dat_chr15_CODAM)
rownames(dat_chr15_CODAM) <- rowns_chr15_CODAM

#LL
dat_chr15_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr15.LL.DS.FORMAT", progress = FALSE)
rowns_chr15_LL <- paste(dat_chr15_LL$CHROM,dat_chr15_LL$POS, sep = ":")
dat_chr15_LL$CHROM <- dat_chr15_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr15_LL),]
colnames(dat_chr15_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr15_LL <- as.matrix(dat_chr15_LL)
rownames(dat_chr15_LL) <- rowns_chr15_LL

#PAN
dat_chr15_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr15.PAN.DS.FORMAT", progress = FALSE)
rowns_chr15_PAN <- paste(dat_chr15_PAN$CHROM, dat_chr15_PAN$POS, sep = ":")
dat_chr15_PAN$CHROM <- dat_chr15_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr15_PAN),]
colnames(dat_chr15_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr15_PAN <- as.matrix(dat_chr15_PAN)
rownames(dat_chr15_PAN) <- rowns_chr15_PAN

#LLS 660Q array
dat_chr15_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr15.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr15_LLS_660Q <- paste(dat_chr15_LLS_660Q$CHROM, dat_chr15_LLS_660Q$POS, sep = ":")
dat_chr15_LLS_660Q$CHROM <- dat_chr15_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr15_LLS_660Q),]
colnames(dat_chr15_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr15_LLS_660Q <- as.matrix(dat_chr15_LLS_660Q)
rownames(dat_chr15_LLS_660Q) <- rowns_chr15_LLS_660Q

#LLS OmniExpr array
dat_chr15_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr15.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr15_LLS_OmniExpr <- paste(dat_chr15_LLS_OmniExpr$CHROM, dat_chr15_LLS_OmniExpr$POS, sep = ":")
dat_chr15_LLS_OmniExpr$CHROM <- dat_chr15_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr15_LLS_OmniExpr),]
colnames(dat_chr15_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr15_LLS_OmniExpr <- as.matrix(dat_chr15_LLS_OmniExpr)
rownames(dat_chr15_LLS_OmniExpr) <- rowns_chr15_LLS_OmniExpr

#RS
dat_chr15_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr15.RS.DS.FORMAT", progress = FALSE)
rowns_chr15_RS <- paste(dat_chr15_RS$CHROM, dat_chr15_RS$POS, sep = ":")
dat_chr15_RS$CHROM <- dat_chr15_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr15_RS),]
colnames(dat_chr15_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr15_RS <- as.matrix(dat_chr15_RS)
rownames(dat_chr15_RS) <- rowns_chr15_RS

#NTR Affy6 array
dat_chr15_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr15.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr15_NTR_Affy6 <- paste(dat_chr15_NTR_Affy6$CHROM,dat_chr15_NTR_Affy6$POS, sep = ":")
dat_chr15_NTR_Affy6$CHROM <- dat_chr15_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr15_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr15_NTR_Affy6),]
colnames(dat_chr15_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr15_NTR_Affy6 <- as.matrix(dat_chr15_NTR_Affy6)
rownames(dat_chr15_NTR_Affy6) <- rowns_chr15_NTR_Affy6

#NTR GoNL array
dat_chr15_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr15.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr15_NTR_GoNL <- paste(dat_chr15_NTR_GoNL$CHROM, dat_chr15_NTR_GoNL$POS, sep = ":")
dat_chr15_NTR_GoNL$CHROM <- dat_chr15_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr15_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr15_NTR_GoNL),]
colnames(dat_chr15_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr15_NTR_GoNL <- as.matrix(dat_chr15_NTR_GoNL)
rownames(dat_chr15_NTR_GoNL) <- rowns_chr15_NTR_GoNL

#keep snps present in every cohort and bind samples
chr15_BIOS_dosages<-list(dat_chr15_CODAM,dat_chr15_LL,dat_chr15_LLS_660Q,dat_chr15_LLS_OmniExpr,dat_chr15_NTR_Affy6,dat_chr15_NTR_GoNL,dat_chr15_PAN,dat_chr15_RS)
rowns <- Reduce(intersect,lapply(chr15_BIOS_dosages, FUN = rownames))
chr15_BIOS_dosages <- lapply(chr15_BIOS_dosages, function(x) x[rowns,])
chr15_BIOS_dosages<-do.call(cbind,chr15_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr15_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr15_BIOS_dosages<-chr15_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr15_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr15_BIOS_dosages <- chr15_BIOS_dosages[, common_columns]
save(chr15_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr15_BIOS_dosages.RData") #yay done

#chr17: dat_chr17_<cohort> -> chr17_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr17_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr17.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr17_CODAM <- paste(dat_chr17_CODAM$CHROM,dat_chr17_CODAM$POS, sep = ":")
dat_chr17_CODAM$CHROM <- dat_chr17_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr17_CODAM),]
colnames(dat_chr17_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr17_CODAM <- as.matrix(dat_chr17_CODAM)
rownames(dat_chr17_CODAM) <- rowns_chr17_CODAM

#LL
dat_chr17_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr17.LL.DS.FORMAT", progress = FALSE)
rowns_chr17_LL <- paste(dat_chr17_LL$CHROM,dat_chr17_LL$POS, sep = ":")
dat_chr17_LL$CHROM <- dat_chr17_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr17_LL),]
colnames(dat_chr17_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr17_LL <- as.matrix(dat_chr17_LL)
rownames(dat_chr17_LL) <- rowns_chr17_LL

#PAN
dat_chr17_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr17.PAN.DS.FORMAT", progress = FALSE)
rowns_chr17_PAN <- paste(dat_chr17_PAN$CHROM, dat_chr17_PAN$POS, sep = ":")
dat_chr17_PAN$CHROM <- dat_chr17_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr17_PAN),]
colnames(dat_chr17_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr17_PAN <- as.matrix(dat_chr17_PAN)
rownames(dat_chr17_PAN) <- rowns_chr17_PAN

#LLS 660Q array
dat_chr17_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr17.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr17_LLS_660Q <- paste(dat_chr17_LLS_660Q$CHROM, dat_chr17_LLS_660Q$POS, sep = ":")
dat_chr17_LLS_660Q$CHROM <- dat_chr17_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr17_LLS_660Q),]
colnames(dat_chr17_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr17_LLS_660Q <- as.matrix(dat_chr17_LLS_660Q)
rownames(dat_chr17_LLS_660Q) <- rowns_chr17_LLS_660Q

#LLS OmniExpr array
dat_chr17_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr17.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr17_LLS_OmniExpr <- paste(dat_chr17_LLS_OmniExpr$CHROM, dat_chr17_LLS_OmniExpr$POS, sep = ":")
dat_chr17_LLS_OmniExpr$CHROM <- dat_chr17_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr17_LLS_OmniExpr),]
colnames(dat_chr17_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr17_LLS_OmniExpr <- as.matrix(dat_chr17_LLS_OmniExpr)
rownames(dat_chr17_LLS_OmniExpr) <- rowns_chr17_LLS_OmniExpr

#RS
dat_chr17_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr17.RS.DS.FORMAT", progress = FALSE)
rowns_chr17_RS <- paste(dat_chr17_RS$CHROM, dat_chr17_RS$POS, sep = ":")
dat_chr17_RS$CHROM <- dat_chr17_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr17_RS),]
colnames(dat_chr17_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr17_RS <- as.matrix(dat_chr17_RS)
rownames(dat_chr17_RS) <- rowns_chr17_RS

#NTR Affy6 array
dat_chr17_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr17.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr17_NTR_Affy6 <- paste(dat_chr17_NTR_Affy6$CHROM,dat_chr17_NTR_Affy6$POS, sep = ":")
dat_chr17_NTR_Affy6$CHROM <- dat_chr17_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr17_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr17_NTR_Affy6),]
colnames(dat_chr17_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr17_NTR_Affy6 <- as.matrix(dat_chr17_NTR_Affy6)
rownames(dat_chr17_NTR_Affy6) <- rowns_chr17_NTR_Affy6

#NTR GoNL array
dat_chr17_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr17.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr17_NTR_GoNL <- paste(dat_chr17_NTR_GoNL$CHROM, dat_chr17_NTR_GoNL$POS, sep = ":")
dat_chr17_NTR_GoNL$CHROM <- dat_chr17_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr17_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr17_NTR_GoNL),]
colnames(dat_chr17_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr17_NTR_GoNL <- as.matrix(dat_chr17_NTR_GoNL)
rownames(dat_chr17_NTR_GoNL) <- rowns_chr17_NTR_GoNL

#keep snps present in every cohort and bind samples
chr17_BIOS_dosages<-list(dat_chr17_CODAM,dat_chr17_LL,dat_chr17_LLS_660Q,dat_chr17_LLS_OmniExpr,dat_chr17_NTR_Affy6,dat_chr17_NTR_GoNL,dat_chr17_PAN,dat_chr17_RS)
rowns <- Reduce(intersect,lapply(chr17_BIOS_dosages, FUN = rownames))
chr17_BIOS_dosages <- lapply(chr17_BIOS_dosages, function(x) x[rowns,])
chr17_BIOS_dosages<-do.call(cbind,chr17_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr17_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr17_BIOS_dosages<-chr17_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr17_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr17_BIOS_dosages <- chr17_BIOS_dosages[, common_columns]
save(chr17_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr17_BIOS_dosages.RData") #yay done

#chr18: dat_chr18_<cohort> -> chr18_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr18_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr18.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr18_CODAM <- paste(dat_chr18_CODAM$CHROM,dat_chr18_CODAM$POS, sep = ":")
dat_chr18_CODAM$CHROM <- dat_chr18_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr18_CODAM),]
colnames(dat_chr18_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr18_CODAM <- as.matrix(dat_chr18_CODAM)
rownames(dat_chr18_CODAM) <- rowns_chr18_CODAM

#LL
dat_chr18_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr18.LL.DS.FORMAT", progress = FALSE)
rowns_chr18_LL <- paste(dat_chr18_LL$CHROM,dat_chr18_LL$POS, sep = ":")
dat_chr18_LL$CHROM <- dat_chr18_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr18_LL),]
colnames(dat_chr18_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr18_LL <- as.matrix(dat_chr18_LL)
rownames(dat_chr18_LL) <- rowns_chr18_LL

#PAN
dat_chr18_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr18.PAN.DS.FORMAT", progress = FALSE)
rowns_chr18_PAN <- paste(dat_chr18_PAN$CHROM, dat_chr18_PAN$POS, sep = ":")
dat_chr18_PAN$CHROM <- dat_chr18_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr18_PAN),]
colnames(dat_chr18_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr18_PAN <- as.matrix(dat_chr18_PAN)
rownames(dat_chr18_PAN) <- rowns_chr18_PAN

#LLS 660Q array
dat_chr18_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr18.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr18_LLS_660Q <- paste(dat_chr18_LLS_660Q$CHROM, dat_chr18_LLS_660Q$POS, sep = ":")
dat_chr18_LLS_660Q$CHROM <- dat_chr18_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr18_LLS_660Q),]
colnames(dat_chr18_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr18_LLS_660Q <- as.matrix(dat_chr18_LLS_660Q)
rownames(dat_chr18_LLS_660Q) <- rowns_chr18_LLS_660Q

#LLS OmniExpr array
dat_chr18_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr18.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr18_LLS_OmniExpr <- paste(dat_chr18_LLS_OmniExpr$CHROM, dat_chr18_LLS_OmniExpr$POS, sep = ":")
dat_chr18_LLS_OmniExpr$CHROM <- dat_chr18_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr18_LLS_OmniExpr),]
colnames(dat_chr18_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr18_LLS_OmniExpr <- as.matrix(dat_chr18_LLS_OmniExpr)
rownames(dat_chr18_LLS_OmniExpr) <- rowns_chr18_LLS_OmniExpr

#RS
dat_chr18_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr18.RS.DS.FORMAT", progress = FALSE)
rowns_chr18_RS <- paste(dat_chr18_RS$CHROM, dat_chr18_RS$POS, sep = ":")
dat_chr18_RS$CHROM <- dat_chr18_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr18_RS),]
colnames(dat_chr18_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr18_RS <- as.matrix(dat_chr18_RS)
rownames(dat_chr18_RS) <- rowns_chr18_RS

#NTR Affy6 array
dat_chr18_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr18.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr18_NTR_Affy6 <- paste(dat_chr18_NTR_Affy6$CHROM,dat_chr18_NTR_Affy6$POS, sep = ":")
dat_chr18_NTR_Affy6$CHROM <- dat_chr18_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr18_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr18_NTR_Affy6),]
colnames(dat_chr18_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr18_NTR_Affy6 <- as.matrix(dat_chr18_NTR_Affy6)
rownames(dat_chr18_NTR_Affy6) <- rowns_chr18_NTR_Affy6

#NTR GoNL array
dat_chr18_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr18.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr18_NTR_GoNL <- paste(dat_chr18_NTR_GoNL$CHROM, dat_chr18_NTR_GoNL$POS, sep = ":")
dat_chr18_NTR_GoNL$CHROM <- dat_chr18_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr18_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr18_NTR_GoNL),]
colnames(dat_chr18_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr18_NTR_GoNL <- as.matrix(dat_chr18_NTR_GoNL)
rownames(dat_chr18_NTR_GoNL) <- rowns_chr18_NTR_GoNL

#keep snps present in every cohort and bind samples
chr18_BIOS_dosages<-list(dat_chr18_CODAM,dat_chr18_LL,dat_chr18_LLS_660Q,dat_chr18_LLS_OmniExpr,dat_chr18_NTR_Affy6,dat_chr18_NTR_GoNL,dat_chr18_PAN,dat_chr18_RS)
rowns <- Reduce(intersect,lapply(chr18_BIOS_dosages, FUN = rownames))
chr18_BIOS_dosages <- lapply(chr18_BIOS_dosages, function(x) x[rowns,])
chr18_BIOS_dosages<-do.call(cbind,chr18_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr18_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr18_BIOS_dosages<-chr18_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr18_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr18_BIOS_dosages <- chr18_BIOS_dosages[, common_columns]
save(chr18_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr18_BIOS_dosages.RData") #yay done

#chr19: dat_chr19_<cohort> -> chr19_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr19_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr19.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr19_CODAM <- paste(dat_chr19_CODAM$CHROM,dat_chr19_CODAM$POS, sep = ":")
dat_chr19_CODAM$CHROM <- dat_chr19_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr19_CODAM),]
colnames(dat_chr19_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr19_CODAM <- as.matrix(dat_chr19_CODAM)
rownames(dat_chr19_CODAM) <- rowns_chr19_CODAM

#LL
dat_chr19_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr19.LL.DS.FORMAT", progress = FALSE)
rowns_chr19_LL <- paste(dat_chr19_LL$CHROM,dat_chr19_LL$POS, sep = ":")
dat_chr19_LL$CHROM <- dat_chr19_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr19_LL),]
colnames(dat_chr19_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr19_LL <- as.matrix(dat_chr19_LL)
rownames(dat_chr19_LL) <- rowns_chr19_LL

#PAN
dat_chr19_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr19.PAN.DS.FORMAT", progress = FALSE)
rowns_chr19_PAN <- paste(dat_chr19_PAN$CHROM, dat_chr19_PAN$POS, sep = ":")
dat_chr19_PAN$CHROM <- dat_chr19_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr19_PAN),]
colnames(dat_chr19_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr19_PAN <- as.matrix(dat_chr19_PAN)
rownames(dat_chr19_PAN) <- rowns_chr19_PAN

#LLS 660Q array
dat_chr19_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr19.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr19_LLS_660Q <- paste(dat_chr19_LLS_660Q$CHROM, dat_chr19_LLS_660Q$POS, sep = ":")
dat_chr19_LLS_660Q$CHROM <- dat_chr19_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr19_LLS_660Q),]
colnames(dat_chr19_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr19_LLS_660Q <- as.matrix(dat_chr19_LLS_660Q)
rownames(dat_chr19_LLS_660Q) <- rowns_chr19_LLS_660Q

#LLS OmniExpr array
dat_chr19_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr19.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr19_LLS_OmniExpr <- paste(dat_chr19_LLS_OmniExpr$CHROM, dat_chr19_LLS_OmniExpr$POS, sep = ":")
dat_chr19_LLS_OmniExpr$CHROM <- dat_chr19_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr19_LLS_OmniExpr),]
colnames(dat_chr19_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr19_LLS_OmniExpr <- as.matrix(dat_chr19_LLS_OmniExpr)
rownames(dat_chr19_LLS_OmniExpr) <- rowns_chr19_LLS_OmniExpr

#RS
dat_chr19_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr19.RS.DS.FORMAT", progress = FALSE)
rowns_chr19_RS <- paste(dat_chr19_RS$CHROM, dat_chr19_RS$POS, sep = ":")
dat_chr19_RS$CHROM <- dat_chr19_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr19_RS),]
colnames(dat_chr19_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr19_RS <- as.matrix(dat_chr19_RS)
rownames(dat_chr19_RS) <- rowns_chr19_RS

#NTR Affy6 array
dat_chr19_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr19.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr19_NTR_Affy6 <- paste(dat_chr19_NTR_Affy6$CHROM,dat_chr19_NTR_Affy6$POS, sep = ":")
dat_chr19_NTR_Affy6$CHROM <- dat_chr19_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr19_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr19_NTR_Affy6),]
colnames(dat_chr19_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr19_NTR_Affy6 <- as.matrix(dat_chr19_NTR_Affy6)
rownames(dat_chr19_NTR_Affy6) <- rowns_chr19_NTR_Affy6

#NTR GoNL array
dat_chr19_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr19.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr19_NTR_GoNL <- paste(dat_chr19_NTR_GoNL$CHROM, dat_chr19_NTR_GoNL$POS, sep = ":")
dat_chr19_NTR_GoNL$CHROM <- dat_chr19_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr19_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr19_NTR_GoNL),]
colnames(dat_chr19_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr19_NTR_GoNL <- as.matrix(dat_chr19_NTR_GoNL)
rownames(dat_chr19_NTR_GoNL) <- rowns_chr19_NTR_GoNL

#keep snps present in every cohort and bind samples
chr19_BIOS_dosages<-list(dat_chr19_CODAM,dat_chr19_LL,dat_chr19_LLS_660Q,dat_chr19_LLS_OmniExpr,dat_chr19_NTR_Affy6,dat_chr19_NTR_GoNL,dat_chr19_PAN,dat_chr19_RS)
rowns <- Reduce(intersect,lapply(chr19_BIOS_dosages, FUN = rownames))
chr19_BIOS_dosages <- lapply(chr19_BIOS_dosages, function(x) x[rowns,])
chr19_BIOS_dosages<-do.call(cbind,chr19_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr19_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr19_BIOS_dosages<-chr19_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr19_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr19_BIOS_dosages <- chr19_BIOS_dosages[, common_columns]
save(chr19_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr19_BIOS_dosages.RData") #yay done

#chr20: dat_chr20_<cohort> -> chr20_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr20_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr20.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr20_CODAM <- paste(dat_chr20_CODAM$CHROM,dat_chr20_CODAM$POS, sep = ":")
dat_chr20_CODAM$CHROM <- dat_chr20_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr20_CODAM),]
colnames(dat_chr20_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr20_CODAM <- as.matrix(dat_chr20_CODAM)
rownames(dat_chr20_CODAM) <- rowns_chr20_CODAM

#LL
dat_chr20_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr20.LL.DS.FORMAT", progress = FALSE)
rowns_chr20_LL <- paste(dat_chr20_LL$CHROM,dat_chr20_LL$POS, sep = ":")
dat_chr20_LL$CHROM <- dat_chr20_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr20_LL),]
colnames(dat_chr20_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr20_LL <- as.matrix(dat_chr20_LL)
rownames(dat_chr20_LL) <- rowns_chr20_LL

#PAN
dat_chr20_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr20.PAN.DS.FORMAT", progress = FALSE)
rowns_chr20_PAN <- paste(dat_chr20_PAN$CHROM, dat_chr20_PAN$POS, sep = ":")
dat_chr20_PAN$CHROM <- dat_chr20_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr20_PAN),]
colnames(dat_chr20_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr20_PAN <- as.matrix(dat_chr20_PAN)
rownames(dat_chr20_PAN) <- rowns_chr20_PAN

#LLS 660Q array
dat_chr20_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr20.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr20_LLS_660Q <- paste(dat_chr20_LLS_660Q$CHROM, dat_chr20_LLS_660Q$POS, sep = ":")
dat_chr20_LLS_660Q$CHROM <- dat_chr20_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr20_LLS_660Q),]
colnames(dat_chr20_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr20_LLS_660Q <- as.matrix(dat_chr20_LLS_660Q)
rownames(dat_chr20_LLS_660Q) <- rowns_chr20_LLS_660Q

#LLS OmniExpr array
dat_chr20_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr20.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr20_LLS_OmniExpr <- paste(dat_chr20_LLS_OmniExpr$CHROM, dat_chr20_LLS_OmniExpr$POS, sep = ":")
dat_chr20_LLS_OmniExpr$CHROM <- dat_chr20_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr20_LLS_OmniExpr),]
colnames(dat_chr20_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr20_LLS_OmniExpr <- as.matrix(dat_chr20_LLS_OmniExpr)
rownames(dat_chr20_LLS_OmniExpr) <- rowns_chr20_LLS_OmniExpr

#RS
dat_chr20_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr20.RS.DS.FORMAT", progress = FALSE)
rowns_chr20_RS <- paste(dat_chr20_RS$CHROM, dat_chr20_RS$POS, sep = ":")
dat_chr20_RS$CHROM <- dat_chr20_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr20_RS),]
colnames(dat_chr20_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr20_RS <- as.matrix(dat_chr20_RS)
rownames(dat_chr20_RS) <- rowns_chr20_RS

#NTR Affy6 array
dat_chr20_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr20.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr20_NTR_Affy6 <- paste(dat_chr20_NTR_Affy6$CHROM,dat_chr20_NTR_Affy6$POS, sep = ":")
dat_chr20_NTR_Affy6$CHROM <- dat_chr20_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr20_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr20_NTR_Affy6),]
colnames(dat_chr20_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr20_NTR_Affy6 <- as.matrix(dat_chr20_NTR_Affy6)
rownames(dat_chr20_NTR_Affy6) <- rowns_chr20_NTR_Affy6

#NTR GoNL array
dat_chr20_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr20.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr20_NTR_GoNL <- paste(dat_chr20_NTR_GoNL$CHROM, dat_chr20_NTR_GoNL$POS, sep = ":")
dat_chr20_NTR_GoNL$CHROM <- dat_chr20_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr20_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr20_NTR_GoNL),]
colnames(dat_chr20_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr20_NTR_GoNL <- as.matrix(dat_chr20_NTR_GoNL)
rownames(dat_chr20_NTR_GoNL) <- rowns_chr20_NTR_GoNL

#keep snps present in every cohort and bind samples
chr20_BIOS_dosages<-list(dat_chr20_CODAM,dat_chr20_LL,dat_chr20_LLS_660Q,dat_chr20_LLS_OmniExpr,dat_chr20_NTR_Affy6,dat_chr20_NTR_GoNL,dat_chr20_PAN,dat_chr20_RS)
rowns <- Reduce(intersect,lapply(chr20_BIOS_dosages, FUN = rownames))
chr20_BIOS_dosages <- lapply(chr20_BIOS_dosages, function(x) x[rowns,])
chr20_BIOS_dosages<-do.call(cbind,chr20_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr20_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr20_BIOS_dosages<-chr20_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr20_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr20_BIOS_dosages <- chr20_BIOS_dosages[, common_columns]
save(chr20_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr20_BIOS_dosages.RData") #yay done

#chr21: dat_chr21_<cohort> -> chr21_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr21_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr21.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr21_CODAM <- paste(dat_chr21_CODAM$CHROM,dat_chr21_CODAM$POS, sep = ":")
dat_chr21_CODAM$CHROM <- dat_chr21_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr21_CODAM),]
colnames(dat_chr21_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr21_CODAM <- as.matrix(dat_chr21_CODAM)
rownames(dat_chr21_CODAM) <- rowns_chr21_CODAM

#LL
dat_chr21_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr21.LL.DS.FORMAT", progress = FALSE)
rowns_chr21_LL <- paste(dat_chr21_LL$CHROM,dat_chr21_LL$POS, sep = ":")
dat_chr21_LL$CHROM <- dat_chr21_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr21_LL),]
colnames(dat_chr21_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr21_LL <- as.matrix(dat_chr21_LL)
rownames(dat_chr21_LL) <- rowns_chr21_LL

#PAN
dat_chr21_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr21.PAN.DS.FORMAT", progress = FALSE)
rowns_chr21_PAN <- paste(dat_chr21_PAN$CHROM, dat_chr21_PAN$POS, sep = ":")
dat_chr21_PAN$CHROM <- dat_chr21_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr21_PAN),]
colnames(dat_chr21_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr21_PAN <- as.matrix(dat_chr21_PAN)
rownames(dat_chr21_PAN) <- rowns_chr21_PAN

#LLS 660Q array
dat_chr21_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr21.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr21_LLS_660Q <- paste(dat_chr21_LLS_660Q$CHROM, dat_chr21_LLS_660Q$POS, sep = ":")
dat_chr21_LLS_660Q$CHROM <- dat_chr21_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr21_LLS_660Q),]
colnames(dat_chr21_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr21_LLS_660Q <- as.matrix(dat_chr21_LLS_660Q)
rownames(dat_chr21_LLS_660Q) <- rowns_chr21_LLS_660Q

#LLS OmniExpr array
dat_chr21_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr21.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr21_LLS_OmniExpr <- paste(dat_chr21_LLS_OmniExpr$CHROM, dat_chr21_LLS_OmniExpr$POS, sep = ":")
dat_chr21_LLS_OmniExpr$CHROM <- dat_chr21_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr21_LLS_OmniExpr),]
colnames(dat_chr21_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr21_LLS_OmniExpr <- as.matrix(dat_chr21_LLS_OmniExpr)
rownames(dat_chr21_LLS_OmniExpr) <- rowns_chr21_LLS_OmniExpr

#RS
dat_chr21_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr21.RS.DS.FORMAT", progress = FALSE)
rowns_chr21_RS <- paste(dat_chr21_RS$CHROM, dat_chr21_RS$POS, sep = ":")
dat_chr21_RS$CHROM <- dat_chr21_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr21_RS),]
colnames(dat_chr21_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr21_RS <- as.matrix(dat_chr21_RS)
rownames(dat_chr21_RS) <- rowns_chr21_RS

#NTR Affy6 array
dat_chr21_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr21.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr21_NTR_Affy6 <- paste(dat_chr21_NTR_Affy6$CHROM,dat_chr21_NTR_Affy6$POS, sep = ":")
dat_chr21_NTR_Affy6$CHROM <- dat_chr21_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr21_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr21_NTR_Affy6),]
colnames(dat_chr21_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr21_NTR_Affy6 <- as.matrix(dat_chr21_NTR_Affy6)
rownames(dat_chr21_NTR_Affy6) <- rowns_chr21_NTR_Affy6

#NTR GoNL array
dat_chr21_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr21.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr21_NTR_GoNL <- paste(dat_chr21_NTR_GoNL$CHROM, dat_chr21_NTR_GoNL$POS, sep = ":")
dat_chr21_NTR_GoNL$CHROM <- dat_chr21_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr21_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr21_NTR_GoNL),]
colnames(dat_chr21_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr21_NTR_GoNL <- as.matrix(dat_chr21_NTR_GoNL)
rownames(dat_chr21_NTR_GoNL) <- rowns_chr21_NTR_GoNL

#keep snps present in every cohort and bind samples
chr21_BIOS_dosages<-list(dat_chr21_CODAM,dat_chr21_LL,dat_chr21_LLS_660Q,dat_chr21_LLS_OmniExpr,dat_chr21_NTR_Affy6,dat_chr21_NTR_GoNL,dat_chr21_PAN,dat_chr21_RS)
rowns <- Reduce(intersect,lapply(chr21_BIOS_dosages, FUN = rownames))
chr21_BIOS_dosages <- lapply(chr21_BIOS_dosages, function(x) x[rowns,])
chr21_BIOS_dosages<-do.call(cbind,chr21_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr21_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr21_BIOS_dosages<-chr21_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr21_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr21_BIOS_dosages <- chr21_BIOS_dosages[, common_columns]
save(chr21_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr21_BIOS_dosages.RData") #yay done

#chr22: dat_chr22_<cohort> -> chr22_BIOS_dosages (common snps, maf > 0.01)
#CODAM
dat_chr22_CODAM <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr22.CODAM.DS.FORMAT", progress = FALSE)
rowns_chr22_CODAM <- paste(dat_chr22_CODAM$CHROM,dat_chr22_CODAM$POS, sep = ":")
dat_chr22_CODAM$CHROM <- dat_chr22_CODAM$POS <- NULL
rownames(mvalues_hrcv1.1_col_CODAM)<-mvalues_hrcv1.1_col_CODAM$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_CODAM<-mvalues_hrcv1.1_col_CODAM[colnames(dat_chr22_CODAM),]
colnames(dat_chr22_CODAM)<-mvalues_hrcv1.1_col_CODAM$uuid
dat_chr22_CODAM <- as.matrix(dat_chr22_CODAM)
rownames(dat_chr22_CODAM) <- rowns_chr22_CODAM

#LL
dat_chr22_LL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr22.LL.DS.FORMAT", progress = FALSE)
rowns_chr22_LL <- paste(dat_chr22_LL$CHROM,dat_chr22_LL$POS, sep = ":")
dat_chr22_LL$CHROM <- dat_chr22_LL$POS <- NULL
rownames(mvalues_hrcv1.1_col_LL)<-mvalues_hrcv1.1_col_LL$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LL<-mvalues_hrcv1.1_col_LL[colnames(dat_chr22_LL),]
colnames(dat_chr22_LL)<-mvalues_hrcv1.1_col_LL$uuid
dat_chr22_LL <- as.matrix(dat_chr22_LL)
rownames(dat_chr22_LL) <- rowns_chr22_LL

#PAN
dat_chr22_PAN <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr22.PAN.DS.FORMAT", progress = FALSE)
rowns_chr22_PAN <- paste(dat_chr22_PAN$CHROM, dat_chr22_PAN$POS, sep = ":")
dat_chr22_PAN$CHROM <- dat_chr22_PAN$POS <- NULL
rownames(mvalues_hrcv1.1_col_PAN)<-mvalues_hrcv1.1_col_PAN$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_PAN<-mvalues_hrcv1.1_col_PAN[colnames(dat_chr22_PAN),]
colnames(dat_chr22_PAN)<-mvalues_hrcv1.1_col_PAN$uuid
dat_chr22_PAN <- as.matrix(dat_chr22_PAN)
rownames(dat_chr22_PAN) <- rowns_chr22_PAN

#LLS 660Q array
dat_chr22_LLS_660Q <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr22.LLS_660Q.DS.FORMAT", progress = FALSE)
rowns_chr22_LLS_660Q <- paste(dat_chr22_LLS_660Q$CHROM, dat_chr22_LLS_660Q$POS, sep = ":")
dat_chr22_LLS_660Q$CHROM <- dat_chr22_LLS_660Q$POS <- NULL
rownames(mvalues_hrcv1.1_col_LLS)<-mvalues_hrcv1.1_col_LLS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_LLS_660Q<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina Human660-Quad")
mvalues_hrcv1.1_col_LLS_660Q<-mvalues_hrcv1.1_col_LLS_660Q[colnames(dat_chr22_LLS_660Q),]
colnames(dat_chr22_LLS_660Q)<-mvalues_hrcv1.1_col_LLS_660Q$uuid
dat_chr22_LLS_660Q <- as.matrix(dat_chr22_LLS_660Q)
rownames(dat_chr22_LLS_660Q) <- rowns_chr22_LLS_660Q

#LLS OmniExpr array
dat_chr22_LLS_OmniExpr <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr22.LLS_OmniExpr.DS.FORMAT", progress = FALSE)
rowns_chr22_LLS_OmniExpr <- paste(dat_chr22_LLS_OmniExpr$CHROM, dat_chr22_LLS_OmniExpr$POS, sep = ":")
dat_chr22_LLS_OmniExpr$CHROM <- dat_chr22_LLS_OmniExpr$POS <- NULL
mvalues_hrcv1.1_col_LLS_OmniExpr<-filter(mvalues_hrcv1.1_col_LLS,gwas_chip=="Illumina OmniExpress")
mvalues_hrcv1.1_col_LLS_OmniExpr<-mvalues_hrcv1.1_col_LLS_OmniExpr[colnames(dat_chr22_LLS_OmniExpr),]
colnames(dat_chr22_LLS_OmniExpr)<-mvalues_hrcv1.1_col_LLS_OmniExpr$uuid
dat_chr22_LLS_OmniExpr <- as.matrix(dat_chr22_LLS_OmniExpr)
rownames(dat_chr22_LLS_OmniExpr) <- rowns_chr22_LLS_OmniExpr

#RS
dat_chr22_RS <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr22.RS.DS.FORMAT", progress = FALSE)
rowns_chr22_RS <- paste(dat_chr22_RS$CHROM, dat_chr22_RS$POS, sep = ":")
dat_chr22_RS$CHROM <- dat_chr22_RS$POS <- NULL
rownames(mvalues_hrcv1.1_col_RS)<-mvalues_hrcv1.1_col_RS$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_RS<-mvalues_hrcv1.1_col_RS[colnames(dat_chr22_RS),]
colnames(dat_chr22_RS)<-mvalues_hrcv1.1_col_RS$uuid
dat_chr22_RS <- as.matrix(dat_chr22_RS)
rownames(dat_chr22_RS) <- rowns_chr22_RS

#NTR Affy6 array
dat_chr22_NTR_Affy6 <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr22.NTR_Affy6.DS.FORMAT", progress = FALSE)
rowns_chr22_NTR_Affy6 <- paste(dat_chr22_NTR_Affy6$CHROM,dat_chr22_NTR_Affy6$POS, sep = ":")
dat_chr22_NTR_Affy6$CHROM <- dat_chr22_NTR_Affy6$POS <- NULL
rownames(mvalues_hrcv1.1_col_NTR)<-mvalues_hrcv1.1_col_NTR$imputation_id_hrcv1.1
mvalues_hrcv1.1_col_NTR_Affy6<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr22_NTR_Affy6))
mvalues_hrcv1.1_col_NTR_Affy6<-mvalues_hrcv1.1_col_NTR_Affy6[colnames(dat_chr22_NTR_Affy6),]
colnames(dat_chr22_NTR_Affy6)<-mvalues_hrcv1.1_col_NTR_Affy6$uuid
dat_chr22_NTR_Affy6 <- as.matrix(dat_chr22_NTR_Affy6)
rownames(dat_chr22_NTR_Affy6) <- rowns_chr22_NTR_Affy6

#NTR GoNL array
dat_chr22_NTR_GoNL <- read_tsv(file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/BIOS_genotype/chr22.NTR_GoNL.DS.FORMAT", progress = FALSE)
rowns_chr22_NTR_GoNL <- paste(dat_chr22_NTR_GoNL$CHROM, dat_chr22_NTR_GoNL$POS, sep = ":")
dat_chr22_NTR_GoNL$CHROM <- dat_chr22_NTR_GoNL$POS <- NULL
mvalues_hrcv1.1_col_NTR_GoNL<-filter(mvalues_hrcv1.1_col_NTR,imputation_id_hrcv1.1%in%colnames(dat_chr22_NTR_GoNL))
mvalues_hrcv1.1_col_NTR_GoNL<-mvalues_hrcv1.1_col_NTR_GoNL[colnames(dat_chr22_NTR_GoNL),]
colnames(dat_chr22_NTR_GoNL)<-mvalues_hrcv1.1_col_NTR_GoNL$uuid
dat_chr22_NTR_GoNL <- as.matrix(dat_chr22_NTR_GoNL)
rownames(dat_chr22_NTR_GoNL) <- rowns_chr22_NTR_GoNL

#keep snps present in every cohort and bind samples
chr22_BIOS_dosages<-list(dat_chr22_CODAM,dat_chr22_LL,dat_chr22_LLS_660Q,dat_chr22_LLS_OmniExpr,dat_chr22_NTR_Affy6,dat_chr22_NTR_GoNL,dat_chr22_PAN,dat_chr22_RS)
rowns <- Reduce(intersect,lapply(chr22_BIOS_dosages, FUN = rownames))
chr22_BIOS_dosages <- lapply(chr22_BIOS_dosages, function(x) x[rowns,])
chr22_BIOS_dosages<-do.call(cbind,chr22_BIOS_dosages)

#drop snps with maf <= 0.01
MAFs <- rowMeans(chr22_BIOS_dosages) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
MAFs <-as.data.frame(MAFs)
MAFs<-filter(MAFs,MAFs>0.01)
chr22_BIOS_dosages<-chr22_BIOS_dosages[rownames(MAFs),]

#keep samples present in mvalues_hrcv1.1
common_columns <- intersect(colnames(chr22_BIOS_dosages), colnames(mvalues_hrcv1.1))
chr22_BIOS_dosages <- chr22_BIOS_dosages[, common_columns]
save(chr22_BIOS_dosages,file = "/exports/archive/kft-pgx-archive/BIOS/HRCv1.1_Imputation/Output_2/chr22_BIOS_dosages.RData") #yay done
