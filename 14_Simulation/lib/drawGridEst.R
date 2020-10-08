################################################################################
# 
### FUNCTION ###################################################################
drawGridEst <- function(est.v=c(), err.v=c(), nboot=1000, SEED=4521){
  
  set.seed(SEED)
  seed.v <- ceiling(runif(n=nboot, min=1, max=1000))
  
  x <- sapply(X=1:nboot, simplify=FALSE, FUN=function(n){
    est.v.sample <- mapply(MEAN=est.v, SD=err.v, FUN=function(MEAN, SD){
      set.seed(seed.v[n])
      rnorm(n=1, mean=MEAN, sd=SD)
    })
    return(est.v.sample)
  })
  
  return(x)
  
}
################################################################################
