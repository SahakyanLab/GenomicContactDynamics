################################################################################
# Make text file of element names in BINREP.MX for A2_HicRepeatExploration.R
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=TRUE)

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
wk.dir = paste0(home.dir, "/DPhil/GCD_polished/18_RepeatVsPersist")

rep.group = "subfam"
binRep.dir = paste0(wk.dir, "/out_RepeatOverlapPerBin/", rep.group)
out.dir = paste0(wk.dir, "/out_makeElementsList/", rep.group)
### OTHER SETTINGS #############################################################
suffix = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chr.v){
  
  load(paste0(binRep.dir, "/", chr, "_BinRep_", suffix, ".RData"))
  element.names <- names(BINREP.MX[1,-(1:3)])
  write(x=element.names, file=paste0(out.dir, "/", suffix, "_", chr, "_elements.txt"),
        ncolumns=1, append=F)
  
}

# rm(list=ls()); gc()