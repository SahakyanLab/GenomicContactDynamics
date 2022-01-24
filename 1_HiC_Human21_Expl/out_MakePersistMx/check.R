################################################################################
# Check if min2Mb PERSIST.MX in persist_HiCNorm_old match with min2Mb PERSIST.MX
# in persist_HiCNorm.
################################################################################
# Set recommended global options

# Avoid left to right partial matching by $
#options(warnPartialMatchDollar=T)

# Expands warnings
#options(warn=1)

wk.dir = "/project/sahakyanlab/ltamon/Database"
old.dir = paste0(wk.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm_old")
new.dir = paste0(wk.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm")
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))

for(chr in chr.v){
  
  print(paste0(chr, "..."))
  
  nme <- paste0(chr, "_Persist_", gcb, ".RData")
  
  load(paste0(old.dir, "/", nme))
  old <- PERSIST.MX
  rm(PERSIST.MX)
  
  load(paste0(new.dir, "/", nme))
  new <- PERSIST.MX
  if(!identical(old, new)){
    warning(paste0(chr, " PERSIST.MX: Not identical"))
  }
  
  rm(PERSIST.MX); gc()
  
}
################################################################################

# rm(list=ls()); gc()