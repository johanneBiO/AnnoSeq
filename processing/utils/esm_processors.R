##########################################################################################
# This script contains functions used to postprocess the attention scores from the ESM-2.
##########################################################################################

library(tidyr)
library(dplyr)
library(jsonlite)
library(stringr)

#---------------------------------------------------------------------------------------#

readJSON <- function(file_path){
  ### Read JSON - each batch is separated by newlines.
  
  data_batch <- readLines(file_path) |>
    lapply(fromJSON)
  
  data_combined <- c()
  
  for (batch in data_batch){
    if(is.list(batch)){
      data_combined <- c(data_combined, batch)
    } else {
      print(as.vector(batch))
      data_combined <- c(data_combined, list(as.vector(batch)))
    }
  }
  
  return(data_combined)
}

readAttnLayer <- function(layer, scr_num, esm_output_path, raw = TRUE){
  ### Read the attention for biological and scrambled sequences for a given layer.
  
  # Define correct file from input.
  layer <- str_pad(layer, width = 2, pad = "0")
  
  if (raw){
    format <- "raw"
  } 
  else {
    format <- "norm"
  }
  
  json_file_bio <- paste(esm_output_path, 
                         "biological_seq/attention/", 
                         format, 
                         "/layer_", 
                         layer, 
                         "_attn_",
                         format,
                         ".json", 
                         sep = "")
  
  # Get the attention from JSON file.
  attn_bio <- readJSON(json_file_bio)
  
  # Get the attention from scrambled sequences.
  attn_scr <- list()
  for (i in 1:scr_num){
    scr <- str_pad(i, width = 2, pad = "0")
    json_file_scr <- paste(esm_output_path, 
                           "scrambled_seq/scramble_", 
                           scr, 
                           "/attention/", 
                           format, 
                           "/layer_", 
                           layer, 
                           "_attn_",
                           format,
                           ".json",
                           sep = "")
    
    attn_scr[[i]] <- readJSON(json_file_scr)
  }
  
  # Transpose list.
  attn_scr <- lapply(1:length(attn_scr[[1]]), function(i) {
    # For each index i (from 1 to 100), extract the i-th element from each list
    lapply(attn_scr, function(inner_list) inner_list[[i]])})
  
  # Save result together.
  attn <- list(bio = attn_bio,
               scr = attn_scr)
  
  return(attn)
}

summarizeScrSignal <- function(scr_list){
  ### Summarize the scrambled signal after aggregating matrices into vectors.
  
  scr_summarized <- lapply(scr_list, function(list) {
    df <- as.data.frame(list)
    colnames(df) <- paste("scr_", 1:dim(df)[2], sep = "")
    df <- df |>
      mutate(position = as.numeric(rownames(df))) 
    df_sum <- df |>
      pivot_longer(cols = -position, 
                   names_to = "scramble", 
                   values_to = "scr_attn") |>
      group_by(position) |>
      summarise(scr_attn_mean = mean(scr_attn),
                scr_attn_sd = sd(scr_attn),
                scr_attn_min = min(scr_attn),
                scr_attn_max = max(scr_attn))
    df <- left_join(df_sum, df, by = "position")
    return(df)})
  
  return(scr_summarized)
}