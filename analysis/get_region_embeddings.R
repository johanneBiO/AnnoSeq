rm(list = ls())

library(here)
library(dplyr)
library(rhdf5)
library(Seurat)

# Read matrix
m_path <- here("data/subset_00100/esm_processed/region_embedding/esm_regions.h5")
mat <- h5read(m_path, "summary_matrix")

# Read metadata
meta_fields <- h5ls(m_path, recursive = TRUE)
meta_cols <- meta_fields$name[meta_fields$group == "/metadata"]

# Construct metadata data.frame
meta_list <- lapply(meta_cols, function(col) {
  h5read(m_path, paste0("metadata/", col))
})

names(meta_list) <- meta_cols
meta_df <- as.data.frame(meta_list, stringsAsFactors = FALSE)

# Ensure rownames match
rownames(mat) <- paste0("dim_", seq_len(nrow(mat)))
colnames(mat) <- paste0("region_", seq_len(ncol(mat)))

rownames(meta_df) <- colnames(mat)

# Build Seurat object
seurat_obj <- CreateSeuratObject(counts = mat)
seurat_obj@meta.data <- meta_df

counts_matrix <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))

seurat_obj <- SetAssayData(seurat_obj, layer = "scale.data", new.data = counts_matrix)
seurat_obj <- RunPCA(seurat_obj, features = rownames(counts_matrix))

ElbowPlot(seurat_obj, reduction = "pca", 30)

seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)

# Plot UMAP
DimPlot(seurat_obj, reduction = "umap", group.by = "feature_type")

#res_path <- here("reports/figures/cluster_embeddings/")
#ggsave(filename = file.path(res_path, "anno_umap_01000.png"), width = 9, height = 5, dpi = 300)
