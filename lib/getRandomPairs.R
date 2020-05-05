################################################################################
# Randomly generate pairs of elements from a set 
# Returns a list of the pairs
# Self pairing not allowed
################################################################################
################################################################################
getRandomPairs <- function(
  pool = bins,
  # Sampling with replacement? see sample()
  replacement = FALSE
){
  pool.len <- 2L
  pair.lst <- list()
  cnt <- 1L
  while( !pool.len%in%c(1,0) ){
    # Tally the pool by identity and sample from names of pool
    # Prevents pairing of bin to itself due to duplicates
    # pool.t <- pool tally
    if(cnt==1L){ pool.t <- table(pool) }
    pair <- sample(x=names(pool.t), size=2, replace=replacement)
    pool.t[ pair[1] ] <- pool.t[ pair[1] ] - 1L
    pool.t[ pair[2] ] <- pool.t[ pair[2] ] - 1L
    pair <- as.numeric(pair)
    if(pair[1] > pair[2]){ pair <- rev(pair) }
    pair.lst[[cnt]] <- pair
    pool.t <- pool.t[pool.t!=0L]
    pool.len <- length(names(pool.t))
    print( pool.len ) 
    cnt <- cnt + 1L
  } # while end
  
  return(pair.lst)
  
}
