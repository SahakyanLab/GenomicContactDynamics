################################################################################
# Count all Hi-C contacts (both short and long-range) per cell line/tissue using 
# MELT.MX.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/2_HiC_Human21_Expl_ext"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/1_Count_contacts"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# MELT.MX directory
melt.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
out.dir = paste0(wk.dir, "/out_count")
### OTHER SETTINGS #############################################################
chr.v <- paste("chr", c(1:22, "X"), sep="") 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
MX <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
  
  # Load MELT.MX
  load(paste0(melt.dir, "/human_", chr, "_allcontacts.RData"))
  
  # Count contacts per cell line/tissue
  out <- apply(X=MELT.MX$upper.tri[,-(1:2)], MARGIN=2, FUN=function(x){
    sum(x!=0)
  })
  
  rm(MELT.MX); gc()
  
  print(paste(chr, " done!"), quote=FALSE)
  
  return(out)
  
})

MX <- do.call("rbind", MX)
MX <- rbind(MX, chrALL=colSums(x=MX))

write.csv(x=MX, file=paste0(out.dir, "/HiCallcontactsCount.csv"),
          row.names=TRUE, quote=FALSE) 

# rm(list=ls())