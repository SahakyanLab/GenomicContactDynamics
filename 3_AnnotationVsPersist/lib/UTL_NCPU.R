#Initialise the NCPU-driven (number of CPUs) setting for do/dopar

################################################################################
###PACKAGES
#CRAN
#install.packages("doMC") 
##package doMC require iterators and parallel packages
#install.packages("iterators") 
#install.packages("parallel") 

################################################################################
if(NCPU > 1){
  #library("iterators")
  #library("parallel")
  library("doMC") # look for the “Parallel” package alternative
  registerDoMC(cores = NCPU)
  `%op%` <- `%dopar%`
  print(paste0("Running with", NCPU, "cores."), quote=F)
} else {
  `%op%` <- `%do%`
}
