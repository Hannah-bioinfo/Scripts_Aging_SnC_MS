# Author: Chuhan LI
# Date: 2025-05-14
# Description: Analysis of transcription factor binding motifs and their regulatory networks

library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

### 1. Prepare DORC Network Components ------------------------
# Calculate k-nearest neighbors for DORC scores
DORC.knn <- FNN::get.knn(
  data = t(scale(Matrix::t(dorcMat.s))), # please refer to 03_Motif_Analysis.R about how to obtain dorcMat.s
  k = dorcK
)$nn.index
rownames(DORC.knn) <- rownames(dorcMat.s)

# Identify DORC-associated peaks
DORCNNpeaks <- unique(dorcTab$Peak[
  dorcTab$Gene %in% c(
    dorcGenes,
    rownames(dorcMat.s)[DORC.knn[rownames(DORC.knn) %in% dorcGenes,]]
  )
])

### 2. Motif Analysis Setup ---------------------------------
# Load human position weight matrices
packagePath <- find.package("FigR", lib.loc = NULL, quiet = TRUE)
pwm <- readRDS(paste0(packagePath, "/data/cisBP_human_pfms_2021.rds"))

# Filter motifs to analyze
motifsToKeep <- intersect(
  names(pwm),
  union(AP1_family, TFsToKeep)
)

# Match motifs to ATAC peaks
motif_position <- motifmatchr::matchMotifs(
  subject = ATAC.se,
  pwms = pwm[motifsToKeep],
  genome = "hg38",
  out = 'positions'
)

### 3. Background Peak Calculation --------------------------
# Prepare ATAC data with GC bias correction
myGenome <- BSgenome.Hsapiens.UCSC.hg38
ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)

# Create background peaks
set.seed(123) # For reproducibility
bg <- chromVAR::getBackgroundPeaks(
  ATAC.se,
  niterations = 50
)

### 4. TSS Annotation Preparation ---------------------------
# Get TSS regions with flanking sequences
TSSg <- FigR::hg38TSSRanges
names(TSSg) <- as.character(TSSg$gene_name)
TSSflank <- GenomicRanges::flank(
  TSSg,
  width = 50000,
  both = TRUE
)

### 5. TF-Gene Regulatory Pair Analysis --------------------
# Example analysis for specific TF-gene pair
tf <- 'JUNB'
sasp <- 'CXCL1'

# Get significant peaks for the TF near the gene
peakSummits <- resize(motif_position[[tf]], width = 1, fix = "center")
TSSflank_SASP <- TSSflank[TSSflank$gene_name == sasp]
overlap_peak_TSS <- overlapsAny(
  query = peakSummits,
  subject = TSSflank_SASP
)

# Create pairs dataframe if overlaps exist
if (length(which(overlap_peak_TSS)) > 0) {
  pairs.df <- as.data.frame(motif_position[[tf]][which(overlap_peak_TSS)])
  
  # Create genomic region string
  region <- paste0(
    pairs.df[1, "seqnames"], 
    ':', 
    min(pairs.df[, "start"]), 
    '-', 
    max(pairs.df[, "end"])
  )
  
  # Visualize the region
  gene_plot <- AnnotationPlot(
    object = MuSC_integrated, 
    region = region
  )
  peak_plot <- PeakPlot(
    object = MuSC_integrated, 
    region = region
  )
  
  # Combine and save plots
  CombineTracks(
    plotlist = list(gene_plot, peak_plot),
    heights = c(1, 1)
  ) + theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 8)
  )
  

### 6. Batch TF Analysis -----------------------------------
# Analyze multiple TFs in batch
AP_Famili <- c('JUN', 'JUNB', 'FOSB', 'FOSL1', 'ATF3', 'BATF', 
             'NFKB1', 'BACH1', 'REL', 'BACH2')
for (tf in AP1_Family) {
  cat("Analyzing TF:", tf, "\n")
  peakSummits <- resize(motif_position[[tf]], width = 1, fix = "center")
  genes <- figR.d.MuSC %>% 
    filter(Motif == tf) %>% 
    filter(Score >= 0.1) %>% 
    pull(DORC) %>% 
    unique()
  if (length(genes) > 0) {
    # Find overlaps with gene regulatory regions
    TSSflank_SASP <- TSSflank[TSSflank$gene_name %in% genes]
    overlap_peak_TSS <- overlapsAny(
      query = peakSummits,
      subject = TSSflank_SASP
    )
  
    if (sum(overlap_peak_TSS) > 0) {
      # Annotate peaks
      pairs.df <- motif_position[[tf]][which(overlap_peak_TSS)]
      peakAnnoList <- annotatePeak(
        pairs.df,
        tssRegion = c(-3000, 3000),
        TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
        flankDistance = 6000
      )
      
      # Save annotation pie chart
      pdf(paste0("./results/", tf, "_peak_annotation.pdf"))
      print(plotAnnoPie(peakAnnoList, cex = 1.5))
      dev.off()
      
    } else {
      cat("No motifs found around genes regulated by", tf, "\n")
    }
  } else {
    cat("No genes significantly regulated by", tf, "\n")
  }
}
