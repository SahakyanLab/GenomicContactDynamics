
################################################################################
# Count all Hi-C contacts (both short and long-range) per tissue per chr using 
# MELT.MX.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon" 
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/2_Expl_contacts")
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/2_Expl_contacts")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(home.dir, "/Database")
# MELT.MX directory
melt.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
out.dir = paste0(wk.dir, "/out_countALLcontacts")
### OTHER SETTINGS #############################################################
chr.v = paste("chr", c(1:22, "X"), sep="") 
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
  
  if( all(rowSums(MELT.MX$upper.tri[,-(1:2)]) > 0) ){
    out <- c(out, allCT=nrow(MELT.MX$upper.tri))
  } else {
    stop(paste0(chr, ": Contact/s not present in any dataset."))
  }
  
  rm(MELT.MX); gc()
  
  print(paste(chr, " done!"), quote=FALSE)
  
  return(out)
  
})

MX <- do.call("rbind", MX)
MX <- rbind(MX, chrALL=colSums(x=MX))

write.csv(x=MX, file=paste0(out.dir, "/HiCallcontactsCount.csv"),
          row.names=TRUE, quote=FALSE) 

# rm(list=ls()); gc()
