library(GSVA)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

# public datasets
CellAge_SnG_Induces = read.table("supplTable2/cellAge_SnG_Induces.txt")[,1]
GenAge_genes = read.table("supplTable2/GenAge_genes.txt")[,1]
SenMayo_genes <- read.table("/lustre/chuhan/1_work/6_human_aging_SASP/00_database/4_SenMayo/SenMayo_genes.txt")[,1]
eigengenes <- read.table("supplTable2/eigengenes_gero.txt")[,1]
sen_sig_genes <- c("CDKN1A", "CDKN2A", "CDKN2B", "CDKN2D", "CDKN1B", "SERPINE1")
geneSets <- list(set1_SenMayo= setdiff(SenMayo_genes, sen_sig_genes), set2_CellAge=setdiff(CellAge_SnG_Induces,sen_sig_genes),
                 set3_GenAge=setdiff(GenAge_genes, sen_sig_genes), set4_eigen = setdiff(eigengenes, sen_sig_genes))
# ssGSVA Score calculation
EC_Mat_orig <- GetAssayData(ECell, slot='data', assay="SCT")
ssGSVA_EC_orig <- gsva(EC_Mat_orig, geneSets, method="ssgsea", ssgsea.norm=TRUE)

ECell$SenMayo_ssGSVA <- ssGSVA_EC_orig['set1_SenMayo',]
ECell$CellAge_ssGSVA <- ssGSVA_EC_orig['set2_CellAge',]
ECell$GenAge_ssGSVA <- ssGSVA_EC_orig['set3_GenAge',]
ECell$eigen_ssGSVA <- ssGSVA_EC_orig['set4_eigen',]
ECell$sen_annotation <- ifelse(ECell$GenAge_ssGSVA > median(ECell$GenAge_ssGSVA)
                                  & ECell$CellAge_ssGSVA > median(ECell$CellAge_ssGSVA)
                                  & ECell$SenMayo_ssGSVA > median(ECell$SenMayo_ssGSVA)
                                  & ECell$eigen_ssGSVA > median(ECell$eigen_ssGSVA),
                                  'SnC', 'non_SnC')
