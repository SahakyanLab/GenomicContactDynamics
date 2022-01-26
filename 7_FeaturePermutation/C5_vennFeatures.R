################################################################################
# Use venn to comment on similaries and differences in enriched/depleted 
# features across regions sets.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
hmap.dir = paste0(wk.dir,"/out_heatmap1/feat_844_raw")
out.dir = paste0(wk.dir,"/out_vennFeatures/feat_844_raw")
### OTHER SETTINGS #############################################################
permtsum.id = "nperm10000_seed662_mxmskfr0"
eval.f.v = c("numOlapA", "comOlap")
hmcol.v = c("Cp21", "CptopCP3", "FC", "ESC", "LC", "hg19")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(venn)
source(paste0(lib, "/doVenn.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(permtsum.id, "_", paste(eval.f.v, collapse="_"))

load(file=paste0(hmap.dir, "/", out.name, "_hmap.RData"))

hvalId.v <- c(`-1`="less", `0`="ns", `1`="greater")
  
for( hval in c(-1,0,1) ){
  
  vennlist <- sapply(X=hmcol.v, simplify=FALSE, FUN=function(x){
    x <- HMAP.MX[,x]
    names(x[x==hval & !is.na(x)])
  })
  
  hval.id <- hvalId.v[as.character(hval)]
  doVenn(vennlist=vennlist, makeVenn=FALSE, saveVenndata=TRUE, venncol="style", 
         filename=paste0(out.dir, "/", hval.id, "_features"))
  
  print(paste0(hval.id, " done!"), quote=FALSE)
  rm(vennlist, hval.id)
  
}

# Find features shared by FC, ESC and LC
HMAP.MX <- HMAP.MX[,c("FC", "ESC", "LC")]
HMAP.MX[!is.na(HMAP.MX)] <- 1
rownames(HMAP.MX)[rowSums(HMAP.MX, na.rm=TRUE)==1]

# rm(list=ls()); gc()



