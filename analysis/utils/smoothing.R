####################################################################################
# This script contains functions used to smooth the attention scores from the ESM-2.
####################################################################################

nearestUneven <- function(x) {
  ### Get the nearest uneven number.
  
  result <- sapply(x, function(val) {
    if (val %% 1 != 0) {
      val <- round(val)
    }
    if (val %% 2 == 0) {
      lower_odd <- val - 1
      upper_odd <- val + 1
      if (abs(val - lower_odd) <= abs(val - upper_odd)) {
        return(lower_odd)
      } else {
        return(upper_odd)
      }
    }
    return(val)
  })
  
  return(result)
}

gaussianKernel <- function(window, sigma){
  ### Generate the Gaussian kernel.
  
  # Generate a symmetric sequence around 0.
  x <- seq(-(window-1)/2, (window-1)/2, by = 1)
  
  # Calculate the Gaussian distribution.
  kernel <- exp(-x^2 / (2 * sigma^2))
  
  # Normalize the kernel so that it sums to 1.
  kernel <- kernel / sum(kernel)
  
  return(kernel)
}

smoothAttention <- function(attention, method = "gaussian", window = 11, sigma = 5){
  ### Smooth the attention signal for one sequence using the mean or gaussian method.
  
  # Check input.
  if (sum(method == c("mean", "gaussian")) != 1){
    stop("Error: Method is not valid. Choose 'mean', 'gaussian' or 'exact'.")
  }
  
  if (window %% 2 == 0){
    stop("Error: Window must be an uneven number.")
  }
  
  length <- length(attention)
  size <- (window-1)/2
  
  # Padding of edges.
  attention_pad <- c(rep(attention[1], size),
                     attention,
                     rep(attention[length], size))
  
  # Apply filter: Moving average.
  if (method == "mean"){
    attention_smooth <- stats::filter(attention_pad,
                                      rep(1/window, window),
                                      sides = 2) |>
      as.numeric()
  }
  
  # Apply filter: Gaussian kernel.
  if (method == "gaussian"){
    kernel <- gaussianKernel(window, sigma)
    attention_smooth <- stats::filter(attention_pad,
                                      kernel,
                                      sides = 2)
  }
  
  # Remove the padding.
  attention_smooth <- attention_smooth[(size+1):(size+length)]
  
  return(attention_smooth)
}
