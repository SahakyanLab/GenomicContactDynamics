################################################################################
# Fraction of long-range contacts per Cp per cell line/tissue per chr
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/2_Expl_contacts"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/1_Count_contacts"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_countLRcontacts")
### OTHER SETTINGS #############################################################
gcb  = "min05Mb"
chr.v = paste0("chr", c(1:22, "X")) 
Cp.v = 1:21
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, "..."), quote=FALSE)
Cp.v <- as.character(sort(Cp.v, decreasing=FALSE))
Cp.v.len <- length(Cp.v)
chr.v.len <- length(chr.v)

for(i in 1:chr.v.len){
  
  chr <- chr.v[i]
  load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  if(i==1){
    ct.v <- setdiff(colnames(PERSIST.MX$hits), c("i", "j"))
    ct.v.len <- length(ct.v)
  }
  mx <- matrix(data=0, nrow=Cp.v.len, ncol=ct.v.len, dimnames=list(Cp.v, ct.v))
  
  for(ct in ct.v){
    ct.TF <- PERSIST.MX$hits[[ct]]>0
    countPCp <- table(PERSIST.MX$ntis[ct.TF])
    countPCp <- countPCp[Cp.v]
    countPCp[is.na(countPCp)] <- 0 # For Cp catergory with no contacts
    mx[names(countPCp),ct] <- countPCp
    rm(countPCp, ct.TF)
  }

  mx <- cbind(mx, allCT=table(PERSIST.MX$ntis)[Cp.v])
  rm(PERSIST.MX); gc()
  mx <- rbind(mx, allCp=colSums(x=mx, na.rm=FALSE))
  write.csv(x=mx, file=paste0(out.dir, "/", gcb, "_", chr, "_HiCLRcontactsCount.csv"),
            row.names=TRUE, quote=FALSE) 
  
  if(i==1){
    LRCOUNT <- mx
  } else {
    LRCOUNT <- LRCOUNT + mx
  }
  rm(mx); gc()
  
  print(paste0(chr, " done!"), quote=FALSE)
  
}

write.csv(x=LRCOUNT, file=paste0(out.dir, "/", gcb, "_chrALL_HiCLRcontactsCount.csv"),
          row.names=TRUE, quote=FALSE) 

# rm(list=ls()); gc()
