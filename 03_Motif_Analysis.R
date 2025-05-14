# DNA Sequence Motif Enrichment and TF Regulation Analysis
# Author: Chuhan LI
# Date: 2025-05-14
# Description: Analysis of transcription factor motifs and their regulatory roles using MuSC as cell type example

# Load required libraries
library(chromVARmotifs)
library(FigR)
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)

### 1. Load Data and Prepare Objects --------------------------
# Load integrated MuSC dataset
MuSC_integrated <- readRDS("./MuSC_integrated.rds")

# Load human PWM motifs
data("human_pwms_v2")

### 2. Motif Analysis Setup -----------------------------------
DefaultAssay(MuSC_integrated) <- 'ATAC'

# Create motif matrix
motif.matrix <- CreateMotifMatrix(
  features = granges(MuSC_integrated),
  pwm = human_pwms_v2,
  genome = 'hg38',
  use.counts = FALSE
)

# Create motif object and add to Seurat object
motif.object <- CreateMotifObject(data = motif.matrix, pwm = human_pwms_v2)
MuSC_integrated <- AddMotifs(
  object = MuSC_integrated,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = human_pwms_v2
)

# Run chromVAR for motif accessibility analysis
MuSC_integrated <- RunChromVAR(
  object = MuSC_integrated,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

### 3. Differential Motif and Gene Analysis ------------------
# Find differentially accessible motifs between senescent and non-senescent cells
DefaultAssay(MuSC_integrated) <- 'chromvar'
motifs_diff <- FindMarkers(
  object = MuSC_integrated,
  ident.1 = "sen",
  ident.2 = "non_sen",
  group.by = 'sen_annotation'
)

# Find differentially expressed genes
DefaultAssay(MuSC_integrated) <- 'SCT'
genes_diff <- FindMarkers(
  object = MuSC_integrated,
  ident.1 = "sen",
  group.by = "sen_annotation",
  only.pos = TRUE,
  logfc.threshold = 0.25
)

### 4. TF-Gene Regulatory Network Analysis -------------------
# Prepare ATAC and RNA data for FigR analysis
ATAC.se <- SummarizedExperiment(
  assays = SimpleList(counts = MuSC_integrated@assays$ATAC@counts),
  colData = MuSC_integrated@meta.data,
  rowRanges = granges(MuSC_integrated)
)

RNAmat <- as.matrix(GetAssayData(MuSC_integrated, assay = "SCT"))

# Run gene-peak correlation analysis
cisCorr <- FigR::runGenePeakcorr(
  ATAC.se = ATAC.se,
  RNAmat = RNAmat,
  genome = "hg38",
  nCores = 15,
  p.cut = 0.05,
  n_bg = 100
)

# Identify DORC genes (Distal Regulatory Elements)
dorcGenes <- dorcJPlot(
  dorcTab = cisCorr,
  cutoff = 1,
  labelTop = 20,
  returnGeneList = TRUE
)

# Calculate DORC scores
dorcMat <- getDORCScores(
  ATAC.se = ATAC.se,
  dorcTab = cisCorr,
  geneList = dorcGenes,
  nCores = 5
)

### 5. TF Regulatory Network Construction -------------------
# Smooth scores using cell KNNs
cellkNN <- FNN::get.knn(
  Embeddings(MuSC_integrated, reduction = "lsi"),
  k = 30
)$nn.index

dorcMat.s <- smoothScoresNN(
  NNmat = cellkNN[,1:30],
  mat = dorcMat,
  nCores = 5
)

# Build TF regulatory network
figr_network <- runFigRGRN(
  ATAC.se = ATAC.se,
  dorcTab = cisCorr,
  genome = "hg38",
  dorcMat = dorcMat.s,
  rnaMat = RNAmat,
  nCores = 5,
  dorcK = 2
)

### 6. Visualization -----------------------------------------
sen_colors <- c('#B1B5BB', '#C220F0') # Non-senescent, senescent
pseudotime_colors <- viridis::viridis(20)

# Plot specific TF motif and gene expression
DefaultAssay(MuSC_integrated) <- 'SCT'
gene_plot <- FeaturePlot(
  MuSC_integrated,
  features = "NFKB1", # taking NFKB1 as an example
  reduction = 'umap.rna',
  cols = sen_colors
)

DefaultAssay(MuSC_integrated) <- 'chromvar'
motif_plot <- FeaturePlot(
  MuSC_integrated,
  features = "NFKB1",
  reduction = 'umap.rna',
  cols = sen_colors
)
