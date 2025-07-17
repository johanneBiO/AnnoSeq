makeSeuratObject <- function(file_path){
  # Read matrix
  embeddings <- h5read(file_path, "summary_matrix")
  
  # Read metadata
  meta_fields <- h5ls(file_path, recursive = TRUE)
  meta_cols <- meta_fields$name[meta_fields$group == "/metadata"]
  
  # Construct metadata data.frame
  meta_list <- lapply(meta_cols, function(col) {
    h5read(file_path, paste0("metadata/", col))
  })
  
  names(meta_list) <- meta_cols
  meta_df <- as.data.frame(meta_list, stringsAsFactors = FALSE)
  
  # Ensure rownames match
  rownames(embeddings) <- paste0("dim_", seq_len(nrow(embeddings)))
  colnames(embeddings) <- paste0("region_", seq_len(ncol(embeddings)))
  rownames(meta_df) <- colnames(embeddings)
  
  # Build Seurat object and include metadata
  seurat_obj <- CreateSeuratObject(counts = embeddings)
  seurat_obj@meta.data <- meta_df
  
  return(seurat_obj)
}

plotHighlight <- function(data, 
                          xvar, 
                          yvar, 
                          highlight_row,
                          x_label = "X", 
                          y_label = "Y",
                          color_line = "#1f618d", 
                          color_point = "#990000",
                          expand_x = 0.05,
                          expand_y = 0.05,
                          n_x,
                          n_y) {
  
  xval <- data[[xvar]]
  yval <- data[[yvar]]
  optimal <- data[highlight_row, , drop = FALSE]
  x_opt <- optimal[[xvar]]
  y_opt <- optimal[[yvar]]
  
  # Compute x/y axis limits with expansion
  xlim <- range(xval)
  x_range <- diff(xlim)
  xlim <- c(xlim[1] - x_range * expand_x, xlim[2] + x_range * expand_x)
  
  ylim <- range(yval)
  y_range <- diff(ylim)
  ylim <- c(ylim[1] - y_range * expand_y, ylim[2] + y_range * expand_y)
  
  p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_path(color = color_line) +
    geom_point(color = color_line, size = 2) +
    
    # Horizontal and vertical guidelines to axis edges
    geom_segment(aes(x = xlim[1], xend = x_opt, y = y_opt, yend = y_opt),
                 color = color_point, linetype = 2) +
    geom_segment(aes(x = x_opt, xend = x_opt, y = ylim[1], yend = y_opt),
                 color = color_point, linetype = 2) +
    
    # Highlighted point
    geom_point(aes(x = x_opt, y = y_opt), color = color_point, size = 5) +
    geom_point(aes(x = x_opt, y = y_opt), color = "#cccccc", size = 2) +
    
    labs(x = x_label, y = y_label) +
    scale_x_continuous(n.breaks = n_x, limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(n.breaks = n_y, limits = ylim, expand = c(0, 0)) +
    theme_minimal() +
    theme(plot.margin = margin(20, 20, 20, 20),
          panel.border = element_rect(color = "black", fill = NA, size = 0.8))
  
  return(p)
}

confusionMat <- function(cluster_res){
  
  # Confusion matrix plot
  conf_mat <- confusionMatrix(cluster_res$cluster_label, cluster_res$feature_type)
  conf_df <- as.data.frame(conf_mat$table) |>
    group_by(Reference) |>
    mutate(Proportion = Freq / sum(Freq)) |>
    ungroup() |>
    mutate(Proportion = case_when(Proportion == 0 ~ NA,
                                  .default = Proportion))
  # Plot the results
  ggplot(conf_df, 
         mapping = aes(x = Prediction, 
                       y = Reference,
                       fill = Proportion)) +
    geom_tile(color = "#999999") +
    geom_text(aes(label = ifelse(Freq > 0, Freq, "")), 
              color = "white",
              size = 3,
              fontface = "bold") +
    scale_fill_gradient(low = "#D2E3EF", 
                        high = "#1f618d",
                        na.value = "#cccccc") +
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
          plot.margin = margin(10, 40, 10, 10))
}

accuracyPerClass <- function(cluster_res, 
                             categories){
  
  set.seed(123)
  
  # Shuffle the true labels
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
  
  # Compute accuracy per class
  accuracy_by_class <- cluster_res |>
    mutate(match = as.character(feature_type) == as.character(cluster_label)) |>
    group_by(feature_type) |>
    summarise(accuracy = sum(match)/length(match)) |>
    left_join(background_accuracy,
              by = "feature_type") |>
    left_join(categories,
              by = "feature_type")
  
  # Plot the results
  ggplot(accuracy_by_class, 
         mapping = aes(x = feature_type,
                       y = accuracy)) +
    geom_bar(aes(y = accuracy), 
             stat = "identity", 
             color = "#1f618d",
             fill = "white",
             alpha = 1,
             size = 1) +
    geom_bar(aes(y = accuracy), 
             stat = "identity", 
             color = "#1f618d",
             fill = "#1f618d",
             alpha = 0.6) +
    geom_bar(aes(y = bg_accuracy), 
             stat = "identity", 
             color = "#1f618d",
             fill = "#990000",
             alpha = 1) +
    geom_text(aes(label = round(accuracy, 2)),
              color = "#1f618d",
              fontface = "bold",
              position = position_dodge(width = 0.8),
              vjust = -0.4,
              size = 3,
              show.legend = FALSE) +
    labs(x = "Annotation Type",
         y = "Accuracy") + 
    facet_grid(cols = vars(category),
               scales = "free_x",
               space = "free_x") +
    scale_y_continuous(limits = c(0, 1), 
                       expand = expansion(mult = c(0, 0), 
                                          add = c(0, 0.05))) +
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
          plot.margin = margin(10, 10, 10, 10))
  
}

adjusted_rand_index <- function(labels_true, 
                                labels_pred) {
  tab <- table(labels_true, labels_pred)
  n <- sum(tab)
  
  # Combinations of counts
  sum_comb_c <- sum(choose(rowSums(tab), 2))
  sum_comb_k <- sum(choose(colSums(tab), 2))
  sum_comb   <- sum(choose(tab, 2))
  
  expected_index <- (sum_comb_c * sum_comb_k) / choose(n, 2)
  max_index <- (sum_comb_c + sum_comb_k) / 2
  ari <- (sum_comb - expected_index) / (max_index - expected_index)
  
  return(ari)
}

normalized_mutual_information <- function(labels_true, 
                                          labels_pred) {
  tab <- table(labels_true, labels_pred)
  joint_prob <- tab / sum(tab)
  
  # Marginal probabilities
  px <- rowSums(joint_prob)
  py <- colSums(joint_prob)
  
  # Mutual Information
  mi <- 0
  for (i in 1:nrow(tab)) {
    for (j in 1:ncol(tab)) {
      if (joint_prob[i, j] > 0) {
        mi <- mi + joint_prob[i, j] * log2(joint_prob[i, j] / (px[i] * py[j]))
      }
    }
  }
  
  # Entropies
  hx <- -sum(px[px > 0] * log2(px[px > 0]))
  hy <- -sum(py[py > 0] * log2(py[py > 0]))
  
  nmi <- mi / sqrt(hx * hy)
  return(nmi)
}
