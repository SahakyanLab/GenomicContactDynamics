################################################################################
# Perform permutation test
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(regioneR)
# library(compiler)
# source(paste0(wk.dir, "/lib/circularRandomizeRegions2.R"))
### FUNCTION ###################################################################
permute <- function(
  out.dir = out.dir, 
  out.name = out.name,
  SEED = SEED,
  # Permutation parameters
  A = A, # Permuted regions. Any format accepted by regioneR::toGRanges.
  B = B, # Any format accepted by regioneR::toGRanges.
  NTIMES= NTIMES,
  genome = genome,
  eval.f.lst = eval.f.lst,
  # Any format accepted by regioneR::toGRanges. No masking if NULL.
  mask.bed = NULL, 
  maxmaskOlapFr = 0,
  nCPU = nCPU,
  # Local z-score calculation parameters. If both NULL, default values are used.
  zs.window = 1e6,
  zs.step = 10000
){
  
  set.seed(SEED)
  PERMT <- permTest(A=A, B=B, ntimes=NTIMES, alternative="auto", 
                    verbose=FALSE, universe=genome, 
                    evaluate.function=eval.f.lst,
                    # Randomisation function
                    randomize.function=regioneR::resampleRegions, 
                    per.chromosome=TRUE, 
                    # Parallel processing; if length(A)<min.parallel - single thread
                    min.parallel=nCPU, 
                    # For reproducibility
                    mc.set.seed=FALSE, mc.cores=nCPU)
  
  if( is.null(zs.window) & is.null(zs.step) ){
    LZSCOR <- localZScore(A=A, B=B, PERMT)
  } else {
    LZSCOR <- localZScore(A=A, B=B, PERMT, window=zs.window, step=zs.step)
  }
  rm(A, B); gc()
  #---------------------------------------Output check
  if( !( identical(names(PERMT), names(LZSCOR)) & identical(names(PERMT), names(eval.f.lst)) ) ){
    stop("PERMT, LZSCOR and eval.f.lst have different orders.")
  }
  
  # Evaluation function that gave zscore=NaN
  ind <- sapply(X=names(LZSCOR), simplify=TRUE, FUN=function(eval.f){
    PERMT[[eval.f]]$zscore
  })
  ind <- which(!is.finite(ind))
  if(length(ind)!=0){
    print(paste0(out.name, ":", length(ind), " NaN/NA zscore/s."), quote=FALSE)
    # Remove NaN/NAs in LZSCOR to avoid error and termination of script
    for(i in ind){
      len <- length( LZSCOR[[i]]$shifted.z.scores )
      # 100 because I just picked a high value
      LZSCOR[[i]]$shifted.z.scores <- rep(100, times=length(LZSCOR[[i]]$shifts))
      print(names(LZSCOR)[i], quote=FALSE)
    }; rm(ind)
  }
  #---------------------------------------Save plots and RData
  pdf(file=paste0(out.dir, "/", out.name, "_plots.pdf"), 
      width=10, height=10)
  plot(PERMT)
  plot(LZSCOR)
  dev.off()
  
  print(paste0("Regions: ", length(A), " A and ", length(B), " B."), quote=FALSE)
  PERMT$region <- c(Anum=length(A), Alen=sum(width(A)), Bnum=length(B), Blen=sum(width(B)))
  save(PERMT, file=paste0(out.dir, "/", out.name, "_permtest.RData"))
  save(LZSCOR, file=paste0(out.dir, "/", out.name, "_zscore.RData"))
  
}
################################################################################
permute <- cmpfun(permute, options=list(suppressUndefined=TRUE))
################################################################################