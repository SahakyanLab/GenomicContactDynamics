################################################################################
# Fraction of unique long-range contact bins per Cp per tissue. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/2_Expl_contacts"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/1_Count_contacts"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_countUbins")
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
    ct.v <- c(setdiff(colnames(PERSIST.MX$hits), c("i", "j")))
    ct.v.len <- length(ct.v)
    
    # Count of all unique contact bins
    allCp <- setNames(obj=rep(NA, times=length(ct.v) + 1),
                      nm=c(ct.v, "allCT"))
  }
  
  mx <- matrix( data=0, nrow=Cp.v.len, ncol=ct.v.len+1L, 
                dimnames=list(Cp.v, c(ct.v, "allCT")) )
  
  for( ct in c(ct.v, "allCT") ){
    if(ct%in%ct.v){
      ct.TF <- PERSIST.MX$hits[[ct]]>0L
    } else if(ct=="allCT"){
      ct.TF <- rep(TRUE, times=length(PERSIST.MX$ntis))
    }
    
    allCp[ct] <- length(unique(unlist(PERSIST.MX$hits[ct.TF,c("i","j")])))
      
    countPCp <- by(data=PERSIST.MX$hits[ct.TF,c("i","j")], INDICES=PERSIST.MX$ntis[ct.TF], 
                   FUN=function(df){
                     return( length(unique(unlist(df))) )
                   })
    countPCp <- countPCp[Cp.v]
    countPCp[is.na(countPCp)] <- 0 # For Cp catergory with no contacts
    mx[names(countPCp),ct] <- countPCp
    rm(countPCp, ct.TF)
  }
  rm(PERSIST.MX); gc()
  mx <- rbind(mx, allCp=allCp)
  write.csv(x=mx, file=paste0(out.dir, "/", gcb, "_", chr, "_countUbins.csv"),
            row.names=TRUE, quote=FALSE) 
  
  if(i==1){
    LRCOUNT <- mx
  } else {
    LRCOUNT <- LRCOUNT + mx
  }
  rm(mx); gc()
  
  print(paste0(chr, " done!"), quote=FALSE)
  
}

write.csv(x=LRCOUNT, file=paste0(out.dir, "/", gcb, "_chrALL_countUbins.csv"),
          row.names=TRUE, quote=FALSE) 

# rm(list=ls()); gc()
