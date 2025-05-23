##########################################################################################
# This script contains functions used to prepare data for classification
##########################################################################################

splitData <- function(data_tbl, 
                      test_size = 0.2, 
                      seed = 10, 
                      test_seq = NULL){
  
  set.seed(seed)
  
  seq <- data_tbl |>
    pull(accession) |>
    unique()
  
  if (is.null(test_seq)){
    test_seq <- sample(seq,
                       size = floor(length(seq) * test_size))
  }
  
  train_seq <- seq[!(seq %in% test_seq)]
  
  data_train <- data_tbl |>
    filter(!(accession %in% test_seq)) |>
    dplyr::select(-accession, -position)
  
  data_test <- data_tbl |>
    filter(accession %in% test_seq) |>
    dplyr::select(-accession, -position)
  
  print(paste("Train positions: ", 
              dim(data_train)[1], 
              " (", 
              round(dim(data_train)[1]/dim(data_tbl)[1], digits = 2), 
              "%)", 
              sep = ""))
  
  print(paste("Test positions: ", 
              dim(data_test)[1], 
              " (", 
              round(dim(data_test)[1]/dim(data_tbl)[1], digits = 2), 
              "%)", 
              sep = ""))
  
  return(list(train = data_train, 
              test = data_test, 
              train_seq = train_seq, 
              test_seq = test_seq))
}

prepForClass <- function(anno_tbl, 
                         attn_tbl, 
                         attn_col, 
                         smooth = TRUE, 
                         w = 11, 
                         s = 2, 
                         wider = TRUE){
  
  # Start time
  start_time <- Sys.time()
  
  attn_col_sym <- sym(attn_col)
  
  # Smoothing in parallel (only if smooth is TRUE)
  if (smooth){
    attn_tbl <- attn_tbl |>
      arrange(accession, position) |>
      group_by(accession, layer) |>
      mutate(!!attn_col_sym := smoothAttention(!!attn_col_sym,
                                               window = w,
                                               sigma = s)) |>
      ungroup()
  }
  
  # Pivot the attention table to wider format
  attn_tbl <- attn_tbl |>
    dplyr::select(accession, position, layer, !!attn_col_sym) |>
    pivot_wider(names_from = layer,
                values_from = !!attn_col_sym,
                names_prefix = "layer_")
  
  # Expand the annotation table.
  anno_tbl <- anno_tbl |>
    arrange(category, feature_type) |>
    distinct(accession, start_position, end_position, category, feature_type) |>
    mutate(position = map2(start_position, end_position, seq)) |>
    unnest(position) |>
    dplyr::select(accession, position, category, feature_type) 
  
  if (isTRUE(wider)){
    # Make annotations into wider format.
    anno_tbl <- anno_tbl |>
      pivot_wider(names_from = feature_type, 
                  values_from = feature_type,
                  values_fn = length,
                  values_fill = 0) |>
      mutate(across(starts_with("feature_"), ~ ifelse(. == 0, 0, 1)),
             across(starts_with("feature_"), as.factor))
    
    # Merge annotation and attention tables, and clean up.
    class_tbl <- attn_tbl |>
      left_join(anno_tbl, by = c("accession", "position")) |>
      mutate(across(everything(), ~ replace_na(.x, 0))) |>
      rowwise() |>
      mutate(feature_uannotated = if_else(sum(c_across(starts_with("feature_"))) == 0, 1, 0)) |>
      rename_with(~ str_remove(.x, "^feature_"), starts_with("feature_"))
    
    # Clean column names (replace spaces and hyphens)
    names(class_tbl) <- gsub("[- ]", "_", names(class_tbl))
  
    } else {
    class_tbl <- attn_tbl |>
      left_join(anno_tbl,
                by = c("accession", "position"),
                relationship = "one-to-many") |>
      mutate(category = case_when(is.na(category) ~ "unannotated",
                                  .default = str_to_lower(category)),
             feature_type = case_when(is.na(feature_type) ~ "unannotated",
                                      .default = str_to_lower(feature_type)),
             category = str_replace_all(category, "[- ]", "_"),
             feature_type = str_replace_all(feature_type, "[- ]", "_"))
  }
  
  # End time
  end_time <- Sys.time()  
  
  elapsed_time <- end_time - start_time
  
  print(paste("Time:", round(elapsed_time, digits = 2)))
  
  return(class_tbl)
}
