# Cell Fate Trajectory Analysis with Monocle3
# Author: Chuhan LI
# Date: 2025-05-4
# Description: Analysis of aging trajectories using pseudotime using MuSC as a cell type example

library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)

### 1. Prepare Data for Monocle --------------------------------
counts.data <- as(as.matrix(MuSC_integrated@assays$RNA@data), 'sparseMatrix')
pheno.data <- new('AnnotatedDataFrame', data = MuSC_integrated@meta.data)
feature.data <- data.frame(gene_short_name = row.names(counts.data), 
                          row.names = row.names(counts.data))
feature.data <- new('AnnotatedDataFrame', data = feature.data)
cds <- newCellDataSet(counts.data,
                     phenoData = pheno.data,
                     featureData = feature.data,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

### 2. Preprocessing and Trajectory Construction --------------
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

# Select ordering genes 
# Example: Using senescence-associated genes
markers_Sen_MuSC <- FindMarkers(object = MuSC_integrated,
                                          ident.1="sen",
                                          group.by="sen_annotation",
                                          features= intersect(union_SnG, rownames(MuSC_integrated)),
                                          logfc.threshold = 0.1)
markers_Sen_MuSC$gene <- rownames(markers_Sen_MuSC)
markers_Sen_MuSC <- arrange(markers_Sen_MuSC, -avg_log2FC)
genes <- markers_Sen_MuSC$gene

# Build trajectory
cds <- setOrderingFilter(cds, genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

# Define root state (youngest cells)
GM_state <- function(cds) {
  if (length(unique(cds$State)) > 1) {
    T0_counts <- table(cds$State, cds$age)[,"Young"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return(1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))

### 3. Visualization ------------------------------------------
sen_colors <- c('#B1B5BB', '#C220F0') # Non-senescent, senescent
pseudotime_colors <- viridis::viridis(20)

# Basic trajectory plot
plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = FALSE) + scale_color_gradientn(colors = pseudotime_colors)

# Grouped by age
plot_cell_trajectory(cds, color_by = "age", show_branch_points = FALSE)

# Grouped by senescence status
plot_cell_trajectory(cds, color_by = "sen_annotation_3", show_branch_points = FALSE) + scale_color_manual(values = sen_colors)

### 4. Pseudotime Analysis ------------------------------------
ordered_pseudot_cds <- cds[, is.finite(cds$Pseudotime)]
a <- as.numeric(ordered_pseudot_cds$Pseudotime)
pseudotime_bin <- cut(a, breaks = seq(min(a), max(a), length.out = 20))
ordered_pseudot_cds$pseudotime_bin <- pseudotime_bin

# Calculate age proportions across pseudotime bins
proportion_df <- data.frame()
for(bin in levels(pseudotime_bin)) {
  meta <- subset(as.data.frame(phenoData(ordered_pseudot_cds)), 
                pseudotime_bin == bin)
  df <- data.frame(table(meta$age)/nrow(meta))
  names(df) <- c("age", "proportion")
  df$bin <- bin
  proportion_df <- rbind(proportion_df, df)
}

# Correlation between pseudotime and aging
proportion_cor <- cor.test(
  subset(proportion_df, age == 'Aged')$bin_num,
  subset(proportion_df, age == 'Aged')$proportion,
  method = 'pearson')

### 5. Differential Expression Analysis -----------------------
# Find markers between trajectory states
## 4 and 5 are two late states identified using Monocle
markers_state <- FindMarkers(object = MuSC_integrated,
                            ident.1 = "4",
                            ident.2 = "5",
                            group.by = "State",
                            logfc.threshold = 0.5)

# Gene-pseudotime correlations
MuSC_gene_pseudotime <- data.frame(gene = rownames(MuSC_Mat_inte), spearman_cor = 0)
for(i in 1:nrow(MuSC_Mat_inte)) {
  res <- cor.test(cds$Pseudotime, MuSC_Mat_inte[i,], 
                 method = "spearman")
  MuSC_gene_pseudotime$spearman_cor[i] <- res$estimate
}

### 6. Functional Enrichment Analysis -------------------------
DEG_late1 <- markers_state %>% filter(avg_log2FC < 0 & p_val_adj < 0.05) %>% pull(gene)
ego_BP_late1 <- enrichGO(gene = DEG_late1,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)

### 7. Heatmap Visualization ----------------------------------
# Subset for late states
MuSC_integrated_late <- subset(MuSC_integrated, subset = (State %in% c('4', '5')))
MuSC_integrated_late$pseudotimec <- ifelse( MuSC_integrated_late$State == '4',
  -MuSC_integrated_late$Pseudotime,
  MuSC_integrated_late$Pseudotime)

# Create heatmap
heatmap_genes <- markers_state %>% filter(abs(avg_log2FC) > 1) %>% pull(gene)
dittoHeatmap(MuSC_integrated_late,
            genes = heatmap_genes,
            annot.by = 'State',
            order.by = 'pseudotimec',
            cluster_cols = FALSE,
            heatmap.colors = colorRampPalette(
              c("#004529", "#f7fcb9", "#f768a1", "#49006a"))(100))
