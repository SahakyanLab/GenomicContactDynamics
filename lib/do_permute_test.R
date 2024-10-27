do_permute_test <- function (control_set, obs_set, olap_set, seed_val = 290, n_iter = 1000) {

  #control_set <- unique(CpGenes$HiC_all)
  #obs_set <- unique(unlist(CpGenes[as.character(1)]))
  #olap_set <- unique(NRgenes)
  
  control_set <- unique(control_set)
  obs_set <- unique(obs_set)
  olap_set <- unique(olap_set)
    
  set.seed(seed_val)  # For reproducibility
  
  # Function to calculate percentage overlap of input with olap_set
  calculate_overlap <- function(input, olap_set) {
    overlap <- length(intersect(input, olap_set))
    return(overlap / length(input) * 100)  # Percentage overlap
  }
  
  observed_overlap <- calculate_overlap(input = obs_set, olap_set = olap_set)
  
  # 2. Permutation test: Randomly permute gene sets and calculate the overlap for each
  permuted_overlaps <- numeric(n_iter)
  obs_set_size <- length(obs_set)
  
  for (i in 1:n_iter) {
    
    permuted_control <- sample(control_set, size = obs_set_size, replace = FALSE)
    permuted_overlaps[i] <- calculate_overlap(input = permuted_control, olap_set = olap_set)
    
  }
  
  # Calculate the p-value for a one-tailed test (observed overlap is lower than most permuted values)
  cat("Permuted Overlap:\n")
  print(summary(permuted_overlaps))
  message("\n")
  message("Observed Overlap:", observed_overlap, "\n")
  
  right_tail_count <- sum(permuted_overlaps >= observed_overlap)
  #if (direction == "lower") {
  if (observed_overlap < mean(permuted_overlaps)) {
    p_value <- sum(permuted_overlaps <= observed_overlap) / n_iter
    names(p_value) <- "one-tailed"
    message("p-value (one-tailed, observed overlap lower):", p_value, "\n")
  #} else if (direction == "higher") {
  } else if (observed_overlap > mean(permuted_overlaps)) {
    p_value <- right_tail_count / n_iter
    names(p_value) <- "one-tailed"
    message("p-value (one-tailed, observed overlap higher):", p_value, "\n")
  } else if (observed_overlap == mean(permuted_overlaps)) {
    left_tail_count <- sum(permuted_overlaps <= -abs(observed_overlap))
    p_value <- (left_tail_count + right_tail_count) / n_iter
    names(p_value) <- "two-tailed"
    message("p-value (two-tailed, observed overlap):", p_value, "\n")
  }
  
  return(p_value)
   
}