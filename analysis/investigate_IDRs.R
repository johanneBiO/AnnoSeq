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
library(org.Hs.eg.db)
library(clusterProfiler)

source("analysis/utils/theme.R")
source("analysis/utils/clust.R")

res_path <- here("reports/figures/cluster_embeddings/")

################################################################################
# DATA IMPORT AND WRANGLING

# Read accessions and annotation
accession <- read_table(here("data/subset_01000/additional/accessions.txt"),  
                        col_names = FALSE)
colnames(accession) <- "accession"

# Read in the Seurat object with true and predicted region embeddings
seurat_comb <- readRDS(file = here("data/subset_01000/esm_processed/region_embedding/seurat_comb.rds"))

# Make region ID as column in metadata
seurat_comb$region_ID <- colnames(seurat_comb)

# Extract the true IDRs and their clusters
seurat_comb$cluster_label <- Idents(seurat_comb)
idr_data <- seurat_comb@meta.data |>
  filter(source == "true",
         feature_type == "Region of interest") |>
  mutate(cluster_number = RNA_snn_res.8) |>
  select(region_ID, accession, category, feature_type, description, cluster_number, cluster_label, start_position, end_position)

# We will only consider clusters also labeled as regions of interest
# We filter out smaller clusters with less than 10 members of potential IDRs
idr_cluster <- idr_data |> 
  filter(cluster_label == "Region of interest") |>
  group_by(cluster_number) |>
  count() |>
  filter(n > 10) |>
  pull(cluster_number)
  
idr_data <- idr_data |>
  filter(cluster_number %in% idr_cluster)

# We write this as a csv file
write_csv(idr_data, here("clustering/results/idr_data.csv"))

################################################################################
# AMINO ACID COMPOSITION

# After computing the amino acid composition, the results are loaded
idr_aa_composition <- read_csv(here("clustering/results/idr_cluster_aa_composition.csv"))
colnames(idr_aa_composition)[1] <- "cluster_number"

idr_aa_composition <- idr_aa_composition |>
  pivot_longer(cols = S:W,
               names_to = "aa",
               values_to = "proportion") |>
  mutate(level2 = case_when(aa %in% c("S", "T", "N", "Q", "C", "Y") ~ "Polar",
                            aa %in% c("K", "R", "H") ~ "Polar (+)",
                            aa %in% c("D", "E") ~ "Polar (-)",
                            aa %in% c("A", "V", "G", "L", "I", "M", "P", "F", "W") ~ "Non-polar"),
         level1 = case_when(level2 == "Non-polar" ~ "Hydrophobic",
                            .default = "Hydrophilic"),
         cluster_number = factor(cluster_number),
         aa = factor(aa),
         level2 = factor(level2, 
                         levels = c("Non-polar", "Polar", "Polar (+)",  "Polar (-)"))) |>
  group_by(cluster_number) |>
  mutate(is_top3 = rank(-proportion, ties.method = "first") <= 3) |>
  ungroup()

# Let's visualize the results
idr_aa_composition |>
ggplot(mapping = aes(x = aa, 
                     y = cluster_number, 
                     fill = proportion)) +
  geom_tile(color = "#999999") +
  geom_text(aes(label = ifelse(is_top3, round(proportion, 2), "")),
            color = "white",
            fontface = "bold",
            size = 3) +
  scale_fill_gradient(low = "white", 
                      high = "#1f618d") +
  facet_grid(cols = vars(level2),
             scales = "free",
             space = "free") +
  labs(x = "Amino Acid",
       y = "Cluster",
       fill = "Proportion") +
  main_theme + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 0, 
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
        plot.margin = margin(10, 10, 10, 10))

ggsave(filename = file.path(res_path, "idr_aa_composition.png"), width = 9, height = 5, dpi = 300)

################################################################################
# COMBINE WITH DISPROT ANNOTATION

# Read disprot results and combine with the IDRs
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

# Let us look at the terms related to the clusters
disprot_go <- disprot |>
  filter(str_detect(term, "GO")) |>
  select(accession, disprot_id, term_namespace, term, term_name, cluster_number) |>
  distinct()

# Let us subset the data and look only at the IDR clusters
subset_idr <- subset(seurat_comb, subset = source == "true")
subset_idr <- subset(subset_idr, comb_info == "Region of interest")

# Define the region identity as the clustering results
Idents(subset_idr) <- subset_idr$RNA_snn_res.8

# Define the 5 clusters to highlight
clusters_to_highlight <- unique(disprot_go$cluster_number) 

# Create a new metadata column to highlight cells that are BOTH in the subset AND belong to those clusters
subset_idr$highlight <- "Other"
cells_to_highlight  <- as.character(subset_idr$RNA_snn_res.8) %in% clusters_to_highlight 
subset_idr$highlight[cells_to_highlight] <- as.character(subset_idr$RNA_snn_res.8[cells_to_highlight])
subset_idr$highlight <- factor(subset_idr$highlight)

# Run the UMAP on the subset
subset_idr <- RunUMAP(subset_idr, 
                      reduction = "pca", 
                      dims = 1:40, 
                      return.model = TRUE)

DimPlot(subset_idr,
        reduction = "umap",
        group.by = "highlight",
        pt.size = 2,
        label = TRUE,
        label.size = 5,
        label.box = TRUE,
        label.color = c("black", "#f39c12", "#7d3c98", "#117864", "#990000", "#1f618d"),
        alpha = 1) +
  scale_color_manual(values = c("#990000", "#1f618d", "#117864", "#f39c12", "#7d3c98", "#cccccc")) +
  scale_fill_manual(values = rep("white", 6)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(filename = file.path(res_path, "idr_umap.png"), width = 7, height = 5, dpi = 300)

# We look at the composition of polar vs. non-polar residues
polar_composition <- idr_aa_composition |>
  group_by(cluster_number, level1) |>
  summarise(p = sum(proportion)) |>
  filter(cluster_number %in% clusters_to_highlight)

################################################################################
# GO ENRICHMENT

# Let us get the gene list for each of the 5 clusters of interest
gene_list <- list()
for (i in 1:5){
  genes <- idr_data |>
    filter(cluster_number == clusters_to_highlight[i]) |>
    pull(accession) |>
    unique()
  
  gene_list[[i]] <- bitr(genes, 
                         fromType = "UNIPROT", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  
  print(paste("Cluster: ", clusters_to_highlight[i], ". Gene members:", length(genes), ". Genes to for GO enrichment: ", dim(gene_list[[i]])[1], ".",
              sep = ""))
}

# We perform the GO enrichment
ego_combined <- data.frame()
for (i in 1:5) {
  ego <- enrichGO(gene = gene_list[[i]]$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  if (!is.null(ego) && nrow(ego) > 0) {
    tmp <- as.data.frame(ego)
    tmp$GeneSet <- paste0("Set ", clusters_to_highlight[i])  # Add gene set label
    ego_combined <- rbind(ego_combined, tmp)
  }
}

# We visualize the results
top_terms <- ego_combined |>
  group_by(GeneSet, ONTOLOGY) |>
  slice_min(order_by = p.adjust, 
            n = 5) |>
  ungroup() |>
  mutate(Description = str_c(str_to_upper(str_sub(Description, 1, 1)),
                             str_sub(Description, 2))) |>
  group_by(GeneSet) |>
  arrange(desc(ONTOLOGY), Count)

top_terms_levels <- unique(top_terms$Description)

top_terms |>
  mutate(Description = factor(Description, 
                              levels = top_terms_levels)) |>
  ggplot(aes(x = Description,
             y = Count,
             fill = p.adjust)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.8),
           alpha = 0.8) +
  scale_fill_gradient(low = "#990000", 
                      high = "#1f618d") +
  facet_grid(rows = vars(GeneSet),
             cols = vars(ONTOLOGY),
             scales = "free",
             space = "free") +
  coord_flip() +
  labs(x = NULL,
       y = "Gene Count",
       fill = "Adj. P-value") +
  scale_y_continuous(breaks = 0:20) +
  main_theme +
  theme(legend.position = "bottom",
        legend.title.align = 0.5,
        strip.background = element_rect(fill="#cccccc", 
                                        color = "#cccccc"),
        strip.text=element_text(color="black", 
                                face = "bold",
                                size = 11),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0.4, "lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(face = "bold", vjust = -1),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 11),
        plot.margin = margin(5, 5, 5, 5)) +
  guides(fill = guide_colorbar(title.position = "bottom", barwidth = 10, barheight = 0.5))

ggsave(filename = file.path(res_path, "idr_goenrich.png"), width = 12, height = 14, dpi = 300)
