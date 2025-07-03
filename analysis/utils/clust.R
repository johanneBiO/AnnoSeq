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

adjusted_rand_index <- function(labels_true, labels_pred) {
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

normalized_mutual_information <- function(labels_true, labels_pred) {
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