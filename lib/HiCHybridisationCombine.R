################################################################################
# Function to combine HYB.MX from all chromosomes
################################################################################
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
### FUNCTION ###################################################################
HiCHybridisationCombine <- function(
  hyb.dir = out.dir,
  persist.dir = persist.dir,
  gcb = "min2Mb",
  chr.v = paste("chr", c(22:1, "X"), sep=""),
  kmer.len = 7,
  affix.hyb = "",
  affix.persist = "",
  saveOut = TRUE
){
  
  patt <- paste0("_Hyb", kmer.len, "_", gcb, affix.hyb, ".RData")
  #path.v <- list.files(path=hyb.dir, pattern=patt, full.names=FALSE)
  
  #chr.v <- gsub(x=path.v, pattern=patt, replacement="")
  #chr.v <- chr.v[chr.v!="chrALL"]
  
  lst <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
    load(file=paste0(hyb.dir, "/", chr, patt))
    print(paste0(chr, " data obtained!"), quote=FALSE)
    return(HYB.MX)
  })
  
  HYB.MX <- do.call("cbind", lst)
  
  # Add Cp row
  ntis.vec <- NULL
  for(chr in chr.v ){
    # Load PERSIST.MX
    load(paste0(persist.dir, "/", chr, "_Persist_", gcb, affix.persist, ".RData"))
    ntis.vec <- c(ntis.vec, PERSIST.MX$ntis)
    rm(PERSIST.MX); gc()
  }
  #rwnme <- dimnames(HYB.MX)[[1]]
  HYB.MX <- rbind(HYB.MX, Cp=ntis.vec)
  #dimnames(HYB.MX)[[1]] <- c(rwnme, "Cp")
    
  rm(lst, ntis.vec, rwnme); gc()
  
  if(saveHYB){
    save(HYB.MX, file=paste0(out.dir, "/chrALL_Hyb", kmer.len, "_", gcb, affix.hyb, ".RData"))
  }
  
  return(HYB.MX)
  
}