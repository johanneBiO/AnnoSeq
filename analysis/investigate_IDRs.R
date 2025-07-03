rm(list = ls())

library(here)
library(dplyr)
library(Seurat)

source("analysis/utils/theme.R")
source("analysis/utils/clust.R")

res_path <- here("reports/figures/cluster_embeddings/")

################################################################################
# DATA IMPORT

# Read accessions and annotation
accession <- read_table(here("data/subset_01000/additional/accessions.txt"),  col_names = FALSE)
colnames(accession) <- "accession"
anno_true <- readRDS("data/complete/annotations/processed/anno.rds") |>
  filter(accession %in% acc$accession)

# Read in the Seurat object with true and predicted region embeddings
seurat_comb <- readRDS(file = here("data/subset_01000/esm_processed/region_embedding/seurat_comb.rds"))

# Extract the true IDRs and their clusters
idr_data <- seurat_comb@meta.data |>
  filter(source == "true",
         description == "Disordered") |>
  mutate(cluster = RNA_snn_res.8) |>
  select(accession, category, feature_type, description, cluster, comb_info, start_position, end_position)

disprot <- read_tsv(file = here("data/external/disprot.tsv"),
                    show_col_types = FALSE) |>
  filter(acc %in% accession$accession) |>
  mutate(accession = acc) |>
  select(-acc) |>
  left_join(idr_data,
            by = c("accession"),
            relationship = "many-to-many") |>
  filter((start == start_position | start > start_position | start >= (start_position - 5)) &
           (end == end_position | end < end_position | end <= (end_position + 5))) |>
  distinct()

test <- as.data.frame(table(disprot$term_name, disprot$cluster)) |>
  filter(Freq != 0, Var1 != "disorder") 

idr_selected <- idr_data |>
  filter(cluster == 21,
         !(accession %in%  disprot$accession))

table(idr_selected$accession)


paste(idr_selected$accession, collapse = " ")

