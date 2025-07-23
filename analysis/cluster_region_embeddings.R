rm(list = ls())

library(here)
library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(rhdf5)
library(caret)
library(Seurat)
library(tidyverse)
library(patchwork)

source("analysis/utils/theme.R")
source("analysis/utils/clust.R")

res_path <- here("reports/figures/cluster_embeddings/")

################################################################################
# DATA IMPORT

# Read accessions and annotation
acc <- read_table(here("data/subset_01000/additional/accessions.txt"),  col_names = FALSE)
colnames(acc) <- "accession"
anno <- readRDS("data/complete/annotations/processed/anno.rds") |>
  filter(accession %in% acc$accession)

# Make Seurat Object
true_path <- here("data/subset_01000/esm_processed/region_embedding/esm_regions_true.h5")
seurat_true <- makeSeuratObject(true_path) 

# Check the dimensions
dim(seurat_true)

# Modify the metadata
feature_type_levels <- anno |>
  select(category, feature_type) |>
  distinct() |>
  arrange(category, feature_type) |>
  pull(feature_type)

feature_type_levels <- paste0(toupper(substring(feature_type_levels, 1, 1)),
                              tolower(substring(feature_type_levels, 2)))
feature_type_levels[feature_type_levels == "Dna-binding region"] <- "DNA-binding region"

seurat_true$feature_type <- as.character(paste0(toupper(substring(seurat_true$feature_type, 1, 1)),
                                                tolower(substring(seurat_true$feature_type, 2))))

seurat_true$feature_type[seurat_true$feature_type == "Dna-binding region"] <- "DNA-binding region"

seurat_true$feature_type <- factor(seurat_true$feature_type,
                                   levels = feature_type_levels)

# Define a category tibble
categories <- anno |>
  group_by(category, feature_type) |>
  count() |>
  select(category, feature_type) |>
  mutate(feature_type = str_c(str_to_upper(str_sub(feature_type, 1, 1)),
                              str_to_lower(str_sub(feature_type, 2))),
         feature_type = case_when(feature_type == "Dna-binding region" ~ "DNA-binding region",
                                  .default = feature_type),
         category = str_c(str_to_upper(str_sub(category, 1, 1)),
                          str_to_lower(str_sub(category, 2))),
         category = case_when(category == "Amino acid modification" ~ "AAM",
                              category == "Molecule processing" ~ "MP",
                              .default = category),
         feature_type = factor(feature_type),
         category = factor(category))

################################################################################
# SCALING

# Avoid normalization
seurat_true <- SetAssayData(seurat_true, 
                            layer = "data", 
                            new.data = GetAssayData(seurat_true, layer = "counts"))

# Scale
all_features <- rownames(seurat_true)
seurat_true <- ScaleData(seurat_true, 
                         features = all_features)

################################################################################
# DIMENSIONALITY REDUCTION

# Run PCA
seurat_true <- RunPCA(seurat_true, 
                      features = all_features, 
                      npcs = 50)

# We evaluate the PCs
stdevs <- seurat_true[["pca"]]@stdev
var_exp <- stdevs^2 / sum(stdevs^2)
cum_var <- cumsum(var_exp)
pc_df <- data.frame(pc = 1:length(cum_var),
                    stdevs = stdevs,
                    cum_var = cum_var,
                    var_exp = var_exp)

# Lets plot the results with 40 PCs marked
p_varexp <- pc_df |>
  filter(pc %% 5 == 0 | pc == 1) |>
  plotHighlight("pc", 
                "cum_var", 
                highlight_row = 9, 
                x_label = "PC", 
                y_label = "Cumulative Explained Variance",
                n_x = 10, 
                n_y = 10) + 
  main_theme 

p_stdev <- pc_df |>
  filter(pc %% 5 == 0 | pc == 1) |>
  plotHighlight("pc", 
                "stdevs", 
                highlight_row = 9, 
                x_label = "PC", 
                y_label = "Standard Deviation",
                n_x = 10, 
                n_y = 10) + 
  main_theme

# We visualize the results in a UMAP plot
seurat_true <- RunUMAP(seurat_true, reduction = "pca", dims = 1:40, return.model = TRUE)

p_true_umap_label <- DimPlot(seurat_true, 
                             reduction = "umap",
                             group.by = "feature_type",
                             cols = palette_19,
                             alpha = 1,
                             pt.size = 1) +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       fill = "Annotation Type") + 
  main_theme +
  theme(plot.title = element_blank())

################################################################################
# CLUSTERING

# Create graph
seurat_true <- FindNeighbors(seurat_true, reduction = "pca", dims = 1:40)

# Clustering:
# The resolution parameter in FindClusters() controls the granularity of clustering.
# Higher resolution = more clusters (finer grouping)
# Lower resolution = fewer clusters (coarser grouping)

# We perform the clustering for a resolution between 0 and 2. 
resolutions <- seq(0.1, 2, by = 0.1)
results <- data.frame(resolution = numeric(), 
                      ARI = numeric(),
                      n_cluster = numeric())

for (res in resolutions) {
  seurat_true <- FindClusters(seurat_true, resolution = res)
  cluster_col <- Idents(seurat_true)
  
  ari <- adjusted_rand_index(seurat_true$feature_type, cluster_col)
  n_cluster <- length(unique(cluster_col))
  
  tab <- table(seurat_true$feature_type, cluster_col)
  
  results <- rbind(results, 
                   data.frame(resolution = res,
                              ARI = ari,
                              n_cluster = n_cluster))
}

# Let's evaluate the results
p_true_ari <- results |>
  filter(resolution < 1.6) |>
  plotHighlight("resolution", 
                "ARI", 
                highlight_row = 5, 
                x_label = "Resolution", 
                y_label = "Adjusted Rand Index",
                n_x = 20, 
                n_y = 10) +
  main_theme

p_true_ncluster <- results |>
  plotHighlight("resolution", 
                "n_cluster", 
                highlight_row = 5, 
                x_label = "Resolution", 
                y_label = "Number of Clusters",
                n_x = 20, 
                n_y = 10) +
  main_theme

################################################################################
# EVALUATE CLUSTERING

tab_df <- as.data.frame(table(cluster = seurat_true$RNA_snn_res.0.5,
                              feature_type = seurat_true$feature_type)) |>
  group_by(feature_type) |>
  mutate(proportion = Freq / sum(Freq)) |>
  ungroup() |>
  left_join(categories,
            by = "feature_type")

p_true_heatmap <- ggplot(tab_df, 
       mapping = aes(x = feature_type,
                     y = cluster,
                     fill = proportion)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "#1f618d") +
  labs(x = "Annotation Type",
       y = "Cluster",
       fill = "Proportion") + 
  facet_grid(cols = vars(category),
             scales = "free_x",
             space = "free_x") +
  main_theme +
  scale_x_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        strip.background = element_rect(fill ="#cccccc", 
                                        color = "#cccccc"),
        strip.text=element_text(color="black", 
                                face = "bold",
                                size = 11),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        plot.margin = margin(20, 20, 20, 20)) 
  
ggsave(filename = file.path(res_path, "true_heatmap_anno_types.png"), width = 6, height = 6, dpi = 300)

# Find the most frequent label per cluster
majority_labels <- seurat_true@meta.data  |>
  group_by(RNA_snn_res.0.5, feature_type) |>
  tally() |>
  top_n(1, wt = n) |>
  group_by(RNA_snn_res.0.5) |>
  ungroup()

# Assign the identity based on the most frequent annotation type within the cluster
Idents(seurat_true) <- seurat_true$RNA_snn_res.0.5
new_cluster_ids <- majority_labels$feature_type
names(new_cluster_ids) <- levels(seurat_true)
seurat_true <- RenameIdents(seurat_true, new_cluster_ids)

# We visualize the results in a UMAP plot
color_index <- which(feature_type_levels %in% unique(Idents(seurat_true)))
p_true_umap_label_new <- DimPlot(seurat_true,
                                 reduction = "umap",
                                 cols = palette_19[color_index],
                                 alpha = 1,
                                 pt.size = 1) +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       fill = "Annotation Type") + 
  main_theme +
  theme(plot.title = element_blank())

# The results are combined
p_true_umap_label <- p_true_umap_label + 
  theme(legend.margin = margin(l = 40))

p_true_umap_label_new <- p_true_umap_label_new +
  theme(legend.position = "none")

((p_varexp + p_true_ari) /
    plot_spacer() / 
    (p_true_umap_label + p_true_umap_label_new)) +
  plot_layout(heights = c(1, 0.1, 2))

ggsave(filename = file.path(res_path, "true_combined_res.png"), width = 14, height = 9, dpi = 300)

# Let's see the proportion of correctly labbeled regions based on the clustering
sum(as.character(Idents(seurat_true)) == as.character(seurat_true$feature_type)) / length(as.character(seurat_true$feature_type))

# Number of represented annotation types after the labeling
unique(as.character(Idents(seurat_true)))

# Percentage of regions belonging to those annotation types not represented after labeling
not_represented <- setdiff(unique(as.character(seurat_true$feature_type)), unique(as.character(Idents(seurat_true)))) 
sum(as.character(seurat_true$feature_type) %in% not_represented) / length(as.character(seurat_true$feature_type))

# Compute the background accuracy
cluster_res <- tibble(cluster = seurat_true$RNA_snn_res.0.5,
                      cluster_label = Idents(seurat_true),
                      feature_type = seurat_true$feature_type)

accuracyPerClass(cluster_res, categories)

ggsave(filename = file.path(res_path, "true_accuracy_per_class.png"), width = 9, height = 5, dpi = 300)

confusionMat(cluster_res)

ggsave(filename = file.path(res_path, "true_confusion_mat.png"), width = 9, height = 7, dpi = 300)

################################################################################
# PCA REPRESENTATION FOR THE PREDICTIONS

# PCA feature loadings (genes × PCs)
loadings <- seurat_true[["pca"]]@feature.loadings

# Load the predicted regions
pred_path <- here("data/subset_01000/esm_processed/region_embedding/esm_regions_pred_1094iter.h5")
seurat_pred <- makeSeuratObject(pred_path)

# Check the dimensions
dim(seurat_pred)

# Scale (without normalization) 
seurat_pred <- SetAssayData(seurat_pred, 
                            layer = "data", 
                            new.data = GetAssayData(seurat_pred, layer = "counts"))

seurat_pred <- ScaleData(seurat_pred, 
                         features = all_features)

# Combine the data
seurat_comb <- seurat_comb <- merge(x = seurat_true,
                                    y = seurat_pred,
                                    add.cell.ids = c("true", "pred"))

# Add the source to metadata
seurat_comb$source <- sub("_.*", "", colnames(seurat_comb))

# Join the layers
seurat_comb[["RNA"]] <- JoinLayers(seurat_comb[["RNA"]])

# Project the data using the loadings from the PCA of the true annotated regions
data_projected <- t(loadings) %*% as.matrix(GetAssayData(seurat_comb, slot = "scale.data")[rownames(loadings), ])

# Add the projected PCA embedding to the PCA slot
seurat_comb[["pca"]] <- CreateDimReducObject(embeddings = t(data_projected),
                                             loadings = seurat_true[["pca"]]@feature.loadings,
                                             key = "PC_",
                                             assay = DefaultAssay(seurat_comb))

################################################################################
# METADATA FOR COMBINED DATA

# Convert Seurat metadata to tibble with cell names as column
# The order is saved as well
meta <- seurat_comb@meta.data |> 
  rownames_to_column("cell_id") |>
  mutate(start_position = as.numeric(start_position),
         end_position = as.numeric(end_position),
         original_order = row_number(),
         source = factor(source, 
                         levels = c("true", "pred"))) |>
  arrange(source, accession, start_position) |>
  mutate(new_order = row_number())

# Make feature type columns comparable with the metadata
anno <- anno |>
  mutate(feature_type = paste0(toupper(substring(feature_type, 1, 1)),
                               tolower(substring(feature_type, 2))),
         feature_type = case_when(feature_type == "Dna-binding region" ~ "DNA-binding region",
                                  .default = feature_type)) |>
  arrange(accession, start_position) |>
  mutate(new_order = row_number())
  
# Join the info
meta_joined <- meta |>
  left_join(anno, 
            by = c("new_order", "accession", "feature_type", "start_position", "end_position")) |>
  select(original_order, cell_id, accession, category, feature_type, description, source, start_position, end_position) |>
  mutate(comb_info = case_when(source == "true" ~ feature_type,
                               source == "pred" ~ "Predicted region"),
         comb_info = factor(comb_info,
                            levels = c("Predicted region", feature_type_levels))) |>
  arrange(original_order) |>
  select(-original_order)

# Reassign metadata (restore rownames)
rownames(meta_joined) <- meta_joined$cell_id
meta_joined$cell_id <- NULL

# Assign back to Seurat object
seurat_comb@meta.data <- meta_joined

################################################################################
# OPTIMIZE CLUSTERING BASED ON DOMAIN SUBTYPES

# We wish to optimize the clustering in order to get the best resolution for 
# identification of annotation subtypes. Domains are very well-studied. 
# Thus, we can optimize based on domains to (hopefully) get the same resolution
# for other annotation types.

# Let's only consider the subtypes with a sufficient number of replicates.
# We choose n >= 5. 
domain_subtypes <- meta_joined |>
  filter(feature_type == "Domain") |>
  group_by(description) |>
  count() |>
  filter(n >= 5)

# Replace NAs - otherwise Seurat removes the regions.
seurat_comb@meta.data[is.na(seurat_comb@meta.data)] <- "unknown"

# Create the graph
seurat_comb <- FindNeighbors(seurat_comb, reduction = "pca", dims = 1:40)

# We perform the clustering for a resolution between 0 and 15. 
resolutions <- seq(1, 15, by = 0.5)
results_domain <- data.frame(resolution = numeric(), 
                             ARI = numeric(), 
                             n_cluster = numeric())

for (res in resolutions) {
  # Perform the clustering on the combined results.
  seurat_comb <- FindClusters(seurat_comb, resolution = res)
  
  # Subset the seurat object and optimize
  seurat_domains <- subset(seurat_comb, subset = description %in% domain_subtypes$description)
  cluster_col <- Idents(seurat_domains)
  
  # Compute metrics
  ari <- adjusted_rand_index(seurat_domains$description, cluster_col)
  ari_acc <- adjusted_rand_index(seurat_domains$accession, cluster_col) 
  n_cluster <- length(unique(cluster_col))
  
  results_domain <- rbind(results_domain, 
                          data.frame(resolution = res, 
                                     ARI = ari,
                                     ARI_ACC = ari_acc,
                                     n_cluster = n_cluster))
}

# We look at the ARI using both domain subtype and sequence accession as reference
results_domain |>
  pivot_longer(cols = ARI:ARI_ACC,
               names_to = "Type",
               values_to = "ARI") |>
  mutate(Reference = case_when(Type == "ARI" ~ "Domain Subtype",
                               Type == "ARI_ACC" ~ "Domain Seq. Acc.")) |>
  ggplot(mapping = aes(x = resolution,
                       y = ARI,
                       color = Reference)) +
  labs(x = "Resolution",
       y = "Adjusted Rand Index") +
  geom_vline(xintercept = 8,
             color = "#cccccc",
             size = 1.5) +
  geom_path() +
  geom_point(size = 2) + 
  scale_color_manual(values = c("#990000", "#1f618d")) + 
  scale_x_continuous(n.breaks = 20, limits = c(1, 15), expand = c(0.05, 0.05)) +
  scale_y_continuous(n.breaks = 8, limits = c(0.05, 0.25), expand = c(0, 0)) +
  main_theme + 
  theme(legend.position = "top",
        plot.margin = margin(5, 5, 5, 10),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8)) +
  guides(fill = guide_legend(ncol = 2,
                             title.position = "top"))

ggsave(filename = file.path(res_path, "ari_domain_subtypes_accession.png"), width = 6, height = 4, dpi = 300)
  
# Evaluate the clustering
tab_df <- as.data.frame(table(cluster = seurat_domains$RNA_snn_res.8,
                              description = seurat_domains$description)) |>
  rename(count = Freq) |>
  group_by(description) |>
  mutate(proportion = count / sum(count),
         proportion = case_when(proportion == 0 ~ NA,
                                .default = proportion)) |>
  ungroup()

domain_clusters <- tab_df |>
  group_by(cluster) |>
  summarise(n = sum(count)) |>
  filter(n > 0)

p_domain <- tab_df |>
  filter(cluster %in% domain_clusters$cluster) |>
  mutate(nonzero = count > 0,
         feature_type = "Domains") |>
  ggplot(aes(x = description,
             y = cluster,
             fill = proportion)) +
  geom_tile(color = "#999999") +
  geom_text(aes(label = ifelse(count > 0, count, "")),
            color = "white",
            fontface = "bold",
            size = 3) +
  scale_fill_gradient2(mid = "#FFD6CC",
                       high = "#990000",
                       na.value = "#cccccc") +
  labs(x = "Domain Subtype",
       y = "Cluster",
       fill = "Proportion") +
  main_theme +
  scale_x_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 50, 
                                   hjust = 1),
        strip.background = element_rect(fill="#cccccc", 
                                        color = "#cccccc"),
        strip.text=element_text(color="black", 
                                face = "bold",
                                size = 11),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.4, "lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(face = "bold", vjust = -1),
        axis.title.y = element_text(face = "bold", vjust = 3),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 11),
        plot.margin = margin(10, 10, 10, 60))

ggsave(filename = file.path(res_path, "true_heatmap_domain_subtypes.png"), width = 9, height = 9, dpi = 300)

################################################################################
# ASSIGN LABELS TO CLUSTERS

# Find the most frequent label per cluster
majority_labels <- seurat_comb@meta.data |>
  filter(source == "true") |>
  group_by(RNA_snn_res.8, feature_type) |>
  tally() |>
  top_n(1, wt = n) |>
  group_by(RNA_snn_res.8) |>
  ungroup()

# Remove cases when there is a tie. 
majority_labels <- majority_labels[-c(11, 61, 106), ]

# We define the identity as the cluster label
Idents(seurat_comb) <- seurat_comb$RNA_snn_res.8
new_cluster_ids <- majority_labels$feature_type
names(new_cluster_ids) <- levels(seurat_comb)
seurat_comb <- RenameIdents(seurat_comb, new_cluster_ids)

# Let's see the proportion of correcly labled regions based on the clustering.
seurat_sub_true <- seurat_comb |>
  subset(subset = source == "true")

sum(as.character(Idents(seurat_sub_true)) == as.character(seurat_sub_true$feature_type)) / length(as.character(seurat_sub_true$feature_type))

# Number of represented annotation types after the labeling
unique(as.character(Idents(seurat_sub_true)))

# Compute the background accuracy
cluster_res_new <- tibble(cluster = seurat_sub_true$RNA_snn_res.8,
                          cluster_label = factor(Idents(seurat_sub_true), 
                                                 levels = feature_type_levels), 
                          feature_type = factor(seurat_sub_true$feature_type,
                                                levels = feature_type_levels))

accuracyPerClass(cluster_res_new, categories)
 
ggsave(filename = file.path(res_path, "true_accuracy_per_class_new.png"), width = 9, height = 5, dpi = 300)

confusionMat(cluster_res_new)

ggsave(filename = file.path(res_path, "true_confusion_mat_new.png"), width = 9, height = 7, dpi = 300)
