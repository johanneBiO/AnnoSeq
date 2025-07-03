rm(list = ls())

library(here)
library(dplyr)
library(rhdf5)
library(Seurat)

source("analysis/utils/theme.R")
source("analysis/utils/clust.R")

res_path <- here("reports/figures/cluster_embeddings/")

################################################################################
# DATA IMPORT

# Read accessions and annotation
acc <- read_table(here("data/subset_01000/additional/accessions.txt"),  col_names = FALSE)
colnames(acc) <- "accession"
anno_true <- readRDS("data/complete/annotations/processed/anno.rds") |>
  filter(accession %in% acc$accession)

# Read in the Seurat object with true and predicted region embeddings
seurat_comb <- readRDS(file = here("data/subset_01000/esm_processed/region_embedding/seurat_comb.rds"))

# Read in accessions and their sequence length to add unannotated sequences
seq_length <- read_csv(file = file.path(here("data/complete/additional/seq_length_sp_cropped.csv"))) |>
  filter(Accession %in% acc$accession)
colnames(seq_length) <- c("accession", "length")

################################################################################
# PREDICTION ACCURACY

# Make a column with the cluster identity
seurat_comb$label <- Idents(seurat_comb)

# Modify the data for the predictions
anno_pred <- seurat_comb@meta.data |>
  filter(source == "pred") |>
  select(accession, feature_type, label, start_position, end_position) |>
  mutate(model = str_remove(feature_type, "model_"),
         model = str_replace_all(model, "_", " "),
         model = str_c(str_to_upper(str_sub(model, 1, 1)),
                       str_to_lower(str_sub(model, 2))),
         feature_type = label) |>
  select(accession, feature_type, model, start_position, end_position)

# Look at the percentile of the detected regions where the model and cluster label overlap.
sum(anno_pred$feature_type == anno_pred$model) / length(anno_pred$model)

# Prepare to add unannotated positions as well
seq_length_unanno <- seq_length |>
  rowwise() |>
  mutate(position = list(1:length)) |>
  unnest(position) |>
  ungroup() |>
  mutate(feature_type = "Unannotated") |>
  select(-length)

# Expand the annotations
anno_true_expanded <- anno_true |>
  mutate(feature_type = str_c(str_to_upper(str_sub(feature_type, 1, 1)),
                              str_to_lower(str_sub(feature_type, 2))),
         position = map2(start_position, end_position, ~ seq(.x, .y))) |>
  unnest(position) |>
  select(accession, position, feature_type)

anno_true_expanded <- bind_rows(anno_true_expanded,
                                anti_join(seq_length_unanno, anno_true_expanded, 
                                          by = c("accession", "position"))) |>
  distinct()

anno_pred_expanded_full <- anno_pred |>
  mutate(position = map2(start_position, end_position, ~ seq(.x, .y))) |>
  unnest(position) |>
  select(accession, position, feature_type, model)

anno_pred_expanded <- anno_pred_expanded_full |>
  select(-model)

anno_pred_expanded <- bind_rows(anno_pred_expanded,
                                anti_join(seq_length_unanno, anno_pred_expanded, 
                                          by = c("accession", "position"))) |>
  distinct()

# Lets look at the number of common rows between the to sets - this is the true positives
correct <- intersect(anno_pred_expanded, anno_true_expanded)

# The number of correctly labled positions (annotated and unannotated) is
dim(correct)[1]/dim(anno_true_expanded)[1]

# Let's find the percentage of annotated positions, that are being detected and correctly annotated
TP <- correct |>
  filter(feature_type != "Unannotated") |>
  nrow()

annotated_pos <- anno_true_expanded |>
  filter(feature_type != "Unannotated") |>
  nrow()

TP/annotated_pos

# Let's evaluate the results only when considering the well-performing models.
good_models <- c("Disulfide bond", "Glycosylation site", "Transmembrane region", "Signal peptide",
                 "Zinc finger region", "Domain", "Region of interest", 
                 "Topological domain", "Active site", "Binding site")

anno_pred_expanded_good_models <- anno_pred_expanded_full |>
  filter(model %in% good_models) |>
  select(-model) |>
  distinct()

anno_true_expanded_good_models <- anno_true_expanded |>
  filter(feature_type %in% good_models)

correct_good_models <- intersect(anno_pred_expanded_good_models, anno_true_expanded_good_models)

# The number of correctly labled positions (annotated and unannotated) is
dim(correct_good_models)[1]/dim(anno_true_expanded_good_models)[1]

################################################################################
# REGION LENGTH DISTRIBUTION

# Compute the region length
anno_true <- anno_true |>
  mutate(region_length = end_position - start_position + 1)

# Compute the region length
anno_pred <- anno_pred |>
  mutate(region_length = end_position - start_position + 1)

print(paste("Average region length (true):", mean(anno_true$region_length)))
print(paste("Average region length (pred):", mean(anno_pred$region_length)))

################################################################################
# PROPORTION OF ANNOTATED POSITIONS IN THE PREDICTIONS

# Modify the data for the predictions
anno_pred <- seurat_comb@meta.data |>
  rownames_to_column(var = "region_id") |>
  filter(source == "pred") |>
  select(region_id, accession, feature_type, label, start_position, end_position, RNA_snn_res.8) |>
  mutate(model = str_remove(feature_type, "model_"),
         model = str_replace_all(model, "_", " "),
         model = str_c(str_to_upper(str_sub(model, 1, 1)),
                       str_to_lower(str_sub(model, 2))),
         feature_type = label,
         position = map2(start_position, end_position, ~ seq(.x, .y))) |>
  unnest(position)

anno_true <- seurat_comb@meta.data |>
  rownames_to_column(var = "region_id") |>
  filter(source == "true") |>
  select(region_id, accession, feature_type, label, start_position, end_position, RNA_snn_res.8) |>
  mutate(model = str_remove(feature_type, "model_"),
         model = str_replace_all(model, "_", " "),
         model = str_c(str_to_upper(str_sub(model, 1, 1)),
                       str_to_lower(str_sub(model, 2))),
         feature_type = label,
         position = map2(start_position, end_position, ~ seq(.x, .y))) |>
  unnest(position)

anno_pred_acc_pos <- anno_pred |>
  select(accession, position)

anno_true_acc_pos <- anno_true |>
  select(accession, position)

common <- intersect(anno_pred_acc_pos, anno_true_acc_pos) |>
  mutate(state = "annotated")

anno_pred <- anno_pred |>
  left_join(common,
            by = c("accession", "position")) |>
  mutate(state = replace_na(state, "unannotated"))

sum(anno_pred$state == "unannotated") + sum(anno_pred$state == "annotated")

proportions <- anno_pred |>
  group_by(region_id) |>
  summarise(total_positions = n(),
            annotated = sum(state == "annotated"),
            proportion_annotated = annotated / total_positions,
            .groups = "drop") |>
  select(region_id, proportion_annotated) |>
  as.data.frame()

rownames(proportions) <- proportions$region_id

seurat_comb <- RunUMAP(seurat_comb, reduction = "pca", dims = 1:40, return.model = TRUE)

subset(seurat_comb, source == "true") |>
  DimPlot(reduction = "umap", 
        cols = palette_19,
        alpha = 0.5,
        pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2", fill = "Annotation Type") + 
  main_theme +
  theme(plot.title = element_blank())

seurat_subset <- subset(seurat_comb, source == "pred")

seurat_subset <- AddMetaData(seurat_subset, metadata = proportions)

seurat_subset |>
  DimPlot(reduction = "umap", 
          group.by = "RNA_snn_res.8",
          label = TRUE,
          #cols = palette_19,
          alpha = 0.5,
          pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2", fill = "Annotation Type") + 
  main_theme +
  theme(plot.title = element_blank())

# Get UMAP coordinates and metadata
umap_data <- Embeddings(seurat_subset, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id")

# Add metadata (make sure it includes `proportion_annotated`)
umap_data$proportion_annotated <- seurat_subset@meta.data$proportion_annotated

# Plot with continuous gradient
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = proportion_annotated)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_gradient(low = "red", high = "lightgray", name = "Proportion Annotated") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  main_theme +
  theme(plot.title = element_blank())

sum(proportions$proportion_annotated == 1)/length(proportions$proportion_annotated)

meta <- seurat_comb@meta.data |>
  mutate(length = end_position-start_position+1)

meta |>
  filter(length < 250) |>
  ggplot(mapping = aes(x = length,
                       colour = source)) +
  geom_density()



quantile(meta$length[meta$source == "pred"])
quantile(meta$length[meta$source == "true"])

meta_region <- seurat_comb@meta.data |>
  filter(description == "Disordered")
   


test <- meta_region |>
  group_by(RNA_snn_res.8, label) |>
  count() |>
  filter(n > 20)


