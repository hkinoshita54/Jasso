# attach libraries ----
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

# Seurat object from matrix, features and barcoodes downloaded from single cell portal ----
expression_matrix <- ReadMtx(
  mtx = "processed_matrix.mtx.gz", features = "processed_features.tsv.gz",
  cells = "processed_barcodes.tsv.gz"
)

meta_data <- read.delim("metaData.txt", stringsAsFactors = F, header = T)
meta_data <- meta_data[-1,] ## remove header
rownames(meta_data) <- meta_data[,1] ## rownames to cell id
meta_data <- meta_data[, c(2,3,7,14:17,21,22)] ## pick up necessary columns

Jasso <- CreateSeuratObject(counts = expression_matrix, meta.data = meta_data)

Jasso[["percent.mt"]] <- PercentageFeatureSet(Jasso, pattern = "^mt-")

# already QC'd
# Jasso <- subset(Jasso, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)



rm(expression_matrix)
rm(meta_data)
gc()

