################################################################################
# Add PREELMTISSDYN.MX from all chromosomes
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = paste0("/Users/ltamon")
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
  } else if (whorunsit == "LiezelCluster"){
    home.dir = paste0("/project/sahakyanlab/ltamon")
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")

rep.group = "subfamALL" # "fam" | "subfam" | "subfam6"
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
PreElmTissDyn.dir = out.dir = paste0(wk.dir, "/out_HicRepeatHeatmapData/", rep.group)
### OTHER SETTINGS #############################################################
# Age rank identifier
out.name = "subfamALL" 
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, "_", rep.group, "..."), quote=FALSE)

col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                    header=TRUE, stringsAsFactors=FALSE)[,col.nme]
agerank <- as.character(agerank)
agerank.len <- length(agerank)

# Initialize final matrix following ELMTISSDYN.MX format
non0cont.MX <- matrix(data=0L, nrow=agerank.len+1L, ncol=21)
dimnames(non0cont.MX) <- list(c("num.contact", agerank), 1:21)

for(chr in chr.v){
  load(file=paste0(PreElmTissDyn.dir, "/", chr, "_", gcb, "_", 
                   out.name, "_PreElmTissDyn.RData"))
  # Check if PreElmTissDyn follows order of agerank; -1 because first column in
  # PreElmTissDyn is ntis
  if( !identical(agerank, rownames(PREELMTISSDYN.MX[-1,])) ){
    stop("Rownames of PREELMTISSDYN.MX not identical to agerank")
  } 
  #for(el in rownames(PREELMTISSDYN.MX)){
  #  non0cont.MX[,el] <- PREELMTISSDYN.MX[,el]
  #}
  non0cont.MX <- non0cont.MX + PREELMTISSDYN.MX
  print(chr)
}
PREELMTISSDYN.MX <- non0cont.MX
rm(non0cont.MX)
# Check if counts of contacts per repeat < total number of contacts
ij.count <- apply(X=PREELMTISSDYN.MX[-1,], MARGIN=1, FUN=function(rw){
  sum(rw>PREELMTISSDYN.MX[1,])
})
ij.count <- sum(ij.count)
if(ij.count!=0){ stop("Checkpoint 2.") }
save(x=PREELMTISSDYN.MX, file=paste0(out.dir, "/chrALL_", gcb, "_", 
                                    out.name, "_PreElmTissDyn.RData"))

# rm(list=ls()); gc()

