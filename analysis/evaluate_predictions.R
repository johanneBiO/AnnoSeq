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
library(VennDiagram)

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

# Define a category tibble
categories <- anno_true |>
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

# Expand the true annotations
anno_true_expanded <- anno_true |>
  mutate(feature_type = str_c(str_to_upper(str_sub(feature_type, 1, 1)),
                              str_to_lower(str_sub(feature_type, 2))),
         position = map2(start_position, end_position, ~ seq(.x, .y))) |>
  unnest(position) |>
  select(accession, position, feature_type)

# Add the remaining unannotated positions
anno_true_expanded <- bind_rows(anno_true_expanded,
                                anti_join(seq_length_unanno, anno_true_expanded, 
                                          by = c("accession", "position"))) |>
  distinct()

# Expand the predictions
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
dim(correct)[1]/dim(anno_pred_expanded)[1]
dim(correct)[1]/dim(anno_true_expanded)[1]

# Let's find the percentage of annotated positions, that are being detected and correctly annotated
TP <- correct |>
  filter(feature_type != "Unannotated") |>
  nrow()

annotated_pred_pos <- anno_pred_expanded |>
  filter(feature_type != "Unannotated") |>
  nrow()

annotated_true_pos <- anno_true_expanded |>
  filter(feature_type != "Unannotated") |>
  nrow()

TP/annotated_pred_pos
TP/annotated_true_pos

# Let's evaluate the results only when considering the well-performing models.
good_models <- c("Disulfide bond", "Glycosylation site", "Transmembrane region", "Signal peptide",
                 "Zinc finger region", "Domain", "Region of interest", 
                 "Topological domain", "Active site", "Binding site")

anno_pred_expanded_good_models <- anno_pred_expanded_full |>
  filter(model %in% good_models) |>
  select(-model) |>
  distinct()

correct_good_models <- intersect(anno_pred_expanded_good_models, anno_true_expanded)

# The number of correctly labbeled positions (annotated and unannotated) is
dim(correct_good_models)[1]/dim(anno_pred_expanded_good_models)[1]

################################################################################
# VENN DIAGRAM

# Lets make a Venn diagram of the two sets of annotated positions
pred_set <- anno_pred_expanded |>
  filter(feature_type != "Unannotated") |>
  mutate(acc_pos = paste(accession, position, sep = "_"),
         acc_pos_anno = paste(accession, position, feature_type, sep = "_")) |>
  distinct()

true_set <- anno_true_expanded |>
  filter(feature_type != "Unannotated") |>
  mutate(acc_pos = paste(accession, position, sep = "_"),
         acc_pos_anno = paste(accession, position, feature_type, sep = "_")) |>
  distinct()

# We first compute the number of distinct positions
length(unique(pred_set$acc_pos))
length(unique(true_set$acc_pos))

# Lets look at how many are in fact annotated among the predicted positions
n_pos_annotated <- intersect(unique(pred_set$acc_pos), unique(true_set$acc_pos))

length(n_pos_annotated)/length(unique(pred_set$acc_pos))

# We make the Venn diagram
venn.diagram(
  x = list(pred_set$acc_pos_anno, true_set$acc_pos_anno),
  category.names = c("Predicted Positions", "True Positions"),
  filename = paste(res_path, "venn.png", sep = ""),
  output=TRUE,
  imagetype="png" ,
  height = 3000 , 
  width = 3000 , 
  resolution = 300,
  compression = "lzw",
  lwd = c(4, 4),
  color = c("#990000", "#1f618d"),
  fill = c("#990000", "#1f618d"),
  alpha = c(0.5, 0.5),
  cex = 2,
  fontface = "bold",
  cat.cex = c(2, 2),
  cat.fontface = c("bold", "bold"))

################################################################################
# INVESTIGATION OF THE SETS

# Find unique and union sets
only_pred <- anti_join(pred_set, true_set)
only_true <- anti_join(true_set, pred_set)
inter <- intersect(pred_set, true_set)

# Add source column
only_pred$source <- "only_pred"
only_true$source <- "only_true"
inter$source <- "inter"

# Combine
set_data <- bind_rows(only_pred, only_true, inter)

# Calculate proportions
proportions <- set_data |>
  group_by(source, feature_type) |>
  summarise(count = n(), 
            .groups = "drop") |>
  group_by(source) |>
  mutate(proportion = count / sum(count),
         feature_type = case_when(feature_type == "Dna-binding region" ~ "DNA-binding region",
                                  .default = feature_type),
         feature_type = factor(feature_type,
                               levels = levels(categories$feature_type)),
         source = factor(source, 
                         levels = c("only_true", "inter", "only_pred")))

# Plot
ggplot(proportions, 
       aes(x = proportion, 
           y = source, 
           fill = feature_type,
           color = feature_type)) +
  labs(title = "",
       fill = "") +
  scale_fill_manual(values = palette_19) + 
  scale_color_manual(values = palette_19) + 
  geom_bar(stat = "identity", 
           position = "fill",
           width = 0.6,
           alpha = 0.5) +
  labs(title = "Class Proportions in Union and Unique Sets",
       y = "Proportion", x = "") +
  theme_void() +
  theme(plot.margin = margin(50, 50, 50, 50)) + 
  guides(color = "none",  # Hide color legend if same as fill
         fill = guide_legend(override.aes = list(alpha = 1))) 
  
ggsave(filename = file.path(res_path, "position_sets_anno_type_proportions.png"), width = 8, height = 4, dpi = 300) 

################################################################################
# POSITION-WISE ACCURACY ON TRUE REGIONS

# Get labels for true regions
anno_true_cluster_label <- seurat_comb@meta.data |>
  filter(source == "true") |>
  select(accession, feature_type, label, start_position, end_position) |>
  mutate(model = str_remove(feature_type, "model_"),
         model = str_replace_all(model, "_", " "),
         model = str_c(str_to_upper(str_sub(model, 1, 1)),
                       str_to_lower(str_sub(model, 2))),
         feature_type = label) |>
  select(accession, feature_type, model, start_position, end_position)

# Expand
anno_true_cluster_label <- anno_true_cluster_label |>
  mutate(position = map2(start_position, end_position, ~ seq(.x, .y))) |>
  unnest(position) |>
  select(accession, position, feature_type) |>
  distinct() |>
  mutate(acc_pos = paste(accession, position, sep = "_"),
         acc_pos_anno = paste(accession, position, feature_type, sep = "_"))

# We first compute the number of distinct positions
length(unique(anno_true_cluster_label$acc_pos))
length(unique(true_set$acc_pos))

# Lets look at how many are in fact annotated among the predicted positions
correct_true <- intersect(unique(anno_true_cluster_label$acc_pos_anno), unique(true_set$acc_pos_anno))

length(correct_true)/length(unique(anno_true_cluster_label$acc_pos_anno))

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
# PROPORTION ANNOTATED POSITIONS IN SEQUENCE

# Proportion annotated in Q8TDD2
anno_true_expanded |>
  filter(accession == "Q8TDD2") |>
  mutate(feature_type = case_when(feature_type != "Unannotated" ~ "Annotated",
                                  .default = feature_type)) |>
  distinct() |>
  group_by(feature_type) |>
  count()

# General
proportion_anno <- anno_true_expanded |>
  mutate(feature_type = case_when(feature_type != "Unannotated" ~ "Annotated",
                                  .default = feature_type)) |>
  distinct() |>
  group_by(accession) |>
  summarise(total_positions = n(),
            anno = sum(feature_type == "Annotated"),
            unanno = sum(feature_type == "Unannotated"),
            proportion = anno / total_positions,
            .groups = "drop") 

proportion_anno |>
  ggplot(mapping = aes(x = proportion)) +
  geom_density(color = "#1f618d",
               fill = "#e0e7ef",
               size = 1) + 
  labs(x = "Proportion Annotated Positions",
       y = "Density") +
  main_theme + 
  scale_x_continuous(limits = c(0, 1), 
                     expand = c(0, 0),
                     labels = function(x) round(x, 2)) + 
  scale_y_continuous(limits = c(0, 1.5),
                     expand = c(0, 0)) + 
  main_theme + 
  theme(plot.margin = margin(20, 20, 20, 20),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8)) 

ggsave(filename = file.path(res_path, "dist_prop_anno.png"), width = 6, height = 4, dpi = 300)

cumulative_unanno <- cumsum(sort(proportion_anno$unanno, decreasing = TRUE))
which(cumulative_unanno >= sum(proportion_anno$unanno) / 2)[1]

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
