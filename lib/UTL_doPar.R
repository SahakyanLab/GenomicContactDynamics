################################################################################
# Initialise the nCPU-driven (number of CPUs) setting for do/dopar
################################################################################
if(nCPU > 1){
  registerDoParallel(cores=nCPU)
  `%op%` <- `%dopar%`
  print(paste0("Running with ", nCPU, " cores."), quote=F)
} else {
  `%op%` <- `%do%`
}