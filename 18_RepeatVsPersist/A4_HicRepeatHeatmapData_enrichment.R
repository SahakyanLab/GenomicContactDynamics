################################################################################
# Make PREELMTISSDYN.MX with bin enrichment per subfamily
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
    fetacp.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/8_FeatureVsPersist/out_FETACP_raw_repeats")
  } else if (whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "subfam"
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
out.dir = paste0(wk.dir, "/out_HicRepeatHeatmapData/subfam_enrichment")
### OTHER SETTINGS #############################################################
fetacp.prefix = "chrALL_min2Mb_hg19"
fetacp.suffix = "repSubfamily_fetacp.RData"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Features
col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                    header=TRUE, stringsAsFactors=FALSE)[,col.nme]; rm(col.nme)
agerank <- as.character(agerank)

agerank.len <- length(agerank)

PREELMTISSDYN.MX <- matrix(data=NA, ncol=21, nrow=agerank.len + 1,
                           dimnames=list(c("num.contact", agerank), as.character(1:21))
                           )

for(ar in agerank){
  
  print(ar, quote=F)
  
  load(paste0(fetacp.dir, "/", fetacp.prefix, "_", ar, "_", fetacp.suffix))
  
  if( identical(dimnames(FETACP.MX)[[2]],  dimnames(FETACP.MX)[[2]]) ){
    PREELMTISSDYN.MX[ar, ] <- unname(FETACP.MX[, "  0"])
  } 

}

PREELMTISSDYN.MX[1,] <- 1

save(PREELMTISSDYN.MX, file=paste0(out.dir, "/chrALL_min2Mb_GiorPubl_PreElmTissDyn.RData"))

# rm(list=ls()); gc()