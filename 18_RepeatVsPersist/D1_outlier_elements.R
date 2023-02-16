################################################################################
# Identify outlying elements based on minrep (shared number)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist"
src.dir = paste0(wk.dir, "/z_ignore_git/out_minRepCounts")
minrep.dir = paste0(src.dir, "/subfamALL_minrep_atleast2sumrep/pairwise")

agerank.file = paste0(wk.dir, "/Repeat_rankingbyAge/repsubfam372.csv")
### OTHER SETTINGS #############################################################
src.id = "chrALL_min2Mb_subfamALL"
minrep.min.median = 1
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
agerank <- read.csv(file=agerank.file, header=T, stringsAsFactors=F)[["repsubfam"]]
agerank <- strsplit(agerank, split=";")[[1]]
agerank.len <- length(agerank)

CHOSEN <- list()
for(ind in 1:agerank.len){
  
  src.suff <- paste0(agerank[[ind]], "_", ind, "_pairwisedifftest")
  
  # minrep
  
  fle <- paste0(minrep.dir, "/", src.id, "_minrepCounts_", src.suff, ".RData")
  if( file.exists(fle) ){
    load(fle)
  } else {
    message(paste0(src.suff, ": skipped." ))
    next
  }
  
  meanmed.mx <- do.call("rbind", TEST$meanmed)
  meanmed.mx[,"MED"]
  
  # Use percentage for proper filtering, look at copynumber and probability of 
  # inserting in a region or contact
  
  
  # 5 transposons
  #CHOSEN[[ind]] <- ifelse(all(meanmed.mx[,"MED"] >= minrep.min.median), agerank[[ind]], NA)
  
  # 8 transposons
  CHOSEN[[ind]] <- ifelse(sum(meanmed.mx[,"MED"] >= minrep.min.median) >= 5, agerank[[ind]], NA)
  
}

CHOSEN <- unlist(CHOSEN)
CHOSEN <- CHOSEN[!is.na(CHOSEN)]
CHOSEN

#> 8 CHOSEN
#[1] "L2"    "MIR3"  "MIRb"  "MIR"   "AluJo" "AluJb" "AluSx" "AluSc"

#> 5 CHOSEN
#[1] "MIR3"  "MIRb"  "MIR"   "AluJb" "AluSx"

# rm(list=ls()); gc()