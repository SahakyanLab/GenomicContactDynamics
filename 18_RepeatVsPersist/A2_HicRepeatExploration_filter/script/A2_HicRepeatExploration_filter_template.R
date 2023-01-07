################################################################################
# Specific for filtering MINREP.MX sumrep based on sumrep value.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" 
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
rep.group = "subfam"
metric = "sumrep"
minelm.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group, "ALL_", metric)
filterMinVal = 2
out.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group, "ALL_", metric, 
                 "_atleast", filterMinVal, "sumrep")
### OTHER SETTINGS #############################################################
suffix = "min2Mb"
chr = "chrCHRREPLACE"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(minelm.dir, "/", chr, "_MinElm_", suffix, ".RData"))
col.len <- length(MINELM.MX[1,])

for(elm.ind in 2:col.len){ # First column is Cp/ntis value
  
  is_failed_filter <- (MINELM.MX[,elm.ind] < filterMinVal) & is.finite(MINELM.MX[,elm.ind])
  MINELM.MX[is_failed_filter, elm.ind] <- NA
  rm(is_failed_filter)
  
}

if( min(MINELM.MX[,-1], na.rm=T) < filterMinVal ){
  stop(paste0(chr, ": Filtering error"))
} else {
  save(MINELM.MX, file=paste0(out.dir, "/", chr, "_MinElm_", suffix, ".RData"))
}

# rm(list=ls()); gc()