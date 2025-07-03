rm(list = ls())

library(here)
library(tidyr)
library(dplyr)
library(readr)
library(rhdf5)
library(Seurat)
library(caret)

source("analysis/utils/theme.R")
source("analysis/utils/clust.R")

res_path <- here("reports/figures/cluster_embeddings/")

palette_19 <- c(
  "#2980b9", "#1f618d", "#B7E4C7", "#1B9E77", "#117864",
  "#DDD6B8", "#CBB67C", "#5C4033",
  "#990000", "#FF0000", "#d35400", "#f39c12", "#F0E442", "black", "#999999", "#cccccc",
  "#E7298A", "#CC79A7", "#7d3c98")

palette_26 <- c(
  "#2980b9", "#1f618d", "blue", "#00CED1", 
  "#B7E4C7", "#1B9E77", "#117864", "green", 
  "#DDD6B8", "#CBB67C", "#8B4513", "#5C4033",
  "#990000", "#FF0000", "#d35400", "#f39c12", "#FFD700" , "#F0E442", 
  "black", "#999999", "#cccccc", 
  "pink", "#E7298A", "#CC79A7", "#7d3c98", "#4B0082")

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

seurat_true$feature_type <- paste0(toupper(substring(seurat_true$feature_type, 1, 1)),
                                   tolower(substring(seurat_true$feature_type, 2)))

seurat_true$feature_type <- factor(seurat_true$feature_type,
                                   levels = feature_type_levels)

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
                y_label = "Variance Explained",
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

p_true_ncluster <-results |>
  plotHighlight("resolution", 
                "n_cluster", 
                highlight_row = 5, 
                x_label = "Resolution", 
                y_label = "Number of Clusters",
                n_x = 20, 
                n_y = 10) +
  main_theme

categories <- anno |>
  group_by(category, feature_type) |>
  count() |>
  select(category, feature_type) |>
  mutate(feature_type = factor(feature_type),
         category = factor(category),
         feature_type = str_c(str_to_upper(str_sub(feature_type, 1, 1)),
                              str_to_lower(str_sub(feature_type, 2))),
         category = str_c(str_to_upper(str_sub(category, 1, 1)),
                          str_to_lower(str_sub(category, 2))),
         category = case_when(category == "Amino acid modification" ~ "AAM",
                              category == "Molecule processing" ~ "MP",
                              .default = category))

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

# Let's see the proportion of correctly labled regions based on the clustering
sum(as.character(Idents(seurat_true)) == as.character(seurat_true$feature_type)) / length(as.character(seurat_true$feature_type))

# Number of represented annotation types after the labeling
unique(as.character(Idents(seurat_true)))

# Compute the background accuracy
cluster_res <- tibble(cluster = seurat_true$RNA_snn_res.0.5,
                      cluster_label = Idents(seurat_true),
                      feature_type = seurat_true$feature_type)

set.seed(123)  

shuffled_accuracies <- map_dfr(1:11, function(i) {
  shuffled_true <- sample(cluster_res$feature_type)
  shuffled_data <- tibble(feature_type = shuffled_true, label = cluster_res$cluster_label)
  
  shuffled_data |>
    mutate(match = as.character(feature_type) == as.character(label)) |>
    group_by(feature_type) |>
    summarise(accuracy = sum(match)/length(match)) |>
    mutate(run = i)
})

# Median background accuracy per class
background_accuracy <- shuffled_accuracies |>
  group_by(feature_type) |>
  summarise(bg_accuracy = median(accuracy))

accuracy_by_class <- cluster_res |>
  mutate(match = as.character(feature_type) == as.character(cluster_label)) |>
  group_by(feature_type) |>
  summarise(accuracy = sum(match)/length(match)) |>
  left_join(background_accuracy,
            by = "feature_type") |>
  left_join(categories,
            by = "feature_type")

ggplot(accuracy_by_class, 
       mapping = aes(x = feature_type)) +
  geom_bar(aes(y = accuracy), 
           stat = "identity", 
           color = "#1f618d",
           fill = "#1f618d",
           alpha = 0.8) +
  geom_bar(aes(y = bg_accuracy), 
           stat = "identity", 
           color = "#1f618d",
           fill = "#990000",
           alpha = 1) +
  labs(x = "Annotation Type",
       y = "Accuracy") + 
  facet_grid(cols = vars(category),
             scales = "free_x",
             space = "free_x") +
  theme_bw() + 
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
        legend.position = "top",
        plot.margin = margin(20, 20, 20, 20))

ggsave(filename = file.path(res_path, "true_accuracy.png"), width = 8, height = 5, dpi = 300)

# Confusion matrix plot
conf_mat <- confusionMatrix(cluster_res$cluster_label, cluster_res$feature_type)
conf_df <- as.data.frame(conf_mat$table) |>
  group_by(Reference) |>
  mutate(Proportion = Freq / sum(Freq)) |>
  ungroup()

ggplot(conf_df, 
       mapping = aes(x = Prediction, 
                     y = Reference, 
                     fill = Proportion)) +
  geom_tile() +
  geom_text(aes(label = Freq), 
            color = "white") +
  scale_fill_gradient(low = "#cccccc", 
                      high = "#1f618d") +
  labs(x = "Predicted Annotation Type",
       y = "True Annotation Type") +
  main_theme + 
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
        plot.margin = margin(20, 20, 20, 20))

ggsave(filename = file.path(res_path, "true_confusion.png"), width = 10, height = 6, dpi = 300)

tibble(start = seurat_true$start_position,
       end = seurat_true$end_position,
       feature_type = seurat_true$feature_type) |>
  mutate(length = as.numeric(end) - as.numeric(start) + 1) |>
  group_by(feature_type) |>
  summarise(mean = mean(length),
            median = median(length)) |>
  left_join(accuracy_by_class,
            by = "feature_type") |>
  ggplot(mapping = aes(x = median,
                       y = accuracy)) + 
  geom_point()

################################################################################
# PCA REPRESENTATION FOR THE PREDICTIONS

# PCA feature loadings (genes × PCs)
loadings <- seurat_true[["pca"]]@feature.loadings

# Load the predicted regions
pred_path <- here("data/subset_01000/esm_processed/region_embedding/esm_regions_pred_1094iter.h5")
seurat_pred <- makeSeuratObject(pred_path)

# Check the dimensionality
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

# Add the projected PCA embeddings to the PCA slot
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
                               tolower(substring(feature_type, 2)))) |>
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

# We perform the clustering for a resolution between 0 and 10. 
resolutions <- seq(1, 10, by = 0.5)
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
  n_cluster <- length(unique(cluster_col))
  
  results_domain <- rbind(results_domain, 
                          data.frame(resolution = res, 
                                     ARI = ari, 
                                     n_cluster = n_cluster))
}

results_domain |>
  plotHighlight("resolution", 
                "ARI", 
                highlight_row = 15, 
                x_label = "Resolution", 
                y_label = "Adjusted Rand Index",
                n_x = 20, 
                n_y = 10) +
  main_theme

results_domain |>
  plotHighlight("resolution", 
                "n_cluster", 
                highlight_row = 15, 
                x_label = "Resolution", 
                y_label = "Number of Domain Clusters",
                n_x = 20, 
                n_y = 10) +
  main_theme

# Evaluate the clustering
tab_df <- as.data.frame(table(cluster = seurat_domains$RNA_snn_res.8,
                              description = seurat_domains$description)) |>
  rename(count = Freq) 

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
             fill = nonzero)) +
  geom_tile(color = "black") +
  geom_text(aes(label = ifelse(count > 0, count, "")),
            color = "white",
            size = 3) +
  scale_fill_manual(
    values = c(`TRUE` = "#1f618d", `FALSE` = "white"),
    labels = c(`TRUE` = "> 0", `FALSE` = "= 0"),
    name = "Count"
  ) +
  labs(x = "Annotation Subtype",
       y = "Cluster") +
  facet_grid(cols = vars(feature_type),
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
        plot.margin = margin(20, 20, 20, 20),
        legend.position = "none") 

ggsave(filename = file.path(res_path, "true_heatmap_domain_subtypes.png"), width = 6, height = 6, dpi = 300)

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

# We look at the UMAP of the true after clustering with the predictions
seurat_comb <- RunUMAP(seurat_comb, reduction = "pca", dims = 1:40, return.model = TRUE)

seurat_comb |>
  subset(subset = source == "true") |>
  DimPlot(reduction = "umap", 
          group.by = "label",
          cols = palette_19,
          alpha = 0.5,
          pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2", fill = "Annotation Type") + 
  main_theme +
  theme(plot.title = element_blank())

# We save the combined data to progress with further analysis in another scripts.
saveRDS(seurat_comb, 
        file = here("data/subset_01000/esm_processed/region_embedding/seurat_comb.rds"))
