#Initialise the nCPU-driven (number of CPUs) setting for do/dopar

################################################################################
################################################################################
#install.packages(doParallel) 

if(nCPU > 1){
  library(doParallel) 
  registerDoParallel(cores=nCPU)
  `%op%` <- `%dopar%`
  print(paste0("Running with", nCPU, "cores."), quote=F)
} else {
  `%op%` <- `%do%`
}
