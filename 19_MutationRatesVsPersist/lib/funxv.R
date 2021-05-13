################################################################################
# Functions for calculating contact values
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(matrixStats)
### FUNCTION ###################################################################
ROWMEANS <- function(x){
  matrixStats::rowMeans2(x=x, na.rm=F)
}
ROWMEDIANS <- function(x){
  matrixStats::rowMedians(x=x, na.rm=F)
}
ROWSDS <- function(x){
  matrixStats::rowSds(x=x, na.rm=F)
}
