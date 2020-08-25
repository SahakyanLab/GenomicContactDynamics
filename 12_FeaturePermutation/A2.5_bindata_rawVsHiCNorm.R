################################################################################
# Check if bin data from HiCNorm values is consistent with bin data from raw 
# values in terms of regions per Cp.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rawbmx.dir = paste0(wk.dir, "/binmx/out_bindata_1perc")
normbmx.dir = paste0(wk.dir, "/binmx/out_bindata_Csmatch_HiCNorm")
#out.dir = paste0(wk.dir, "/binmx/out_bindata_Csmatch_HiCNorm")
### OTHER SETTINGS #############################################################
gcb = "min2Mb" # "min2Mb" | "min05Mb"
chr.v = paste0("chr", c(1:22, "X"), sep="")
Cs.cutoff.rw = 0.01
Cs.cutoff.norm = 2
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0("Cs cutoff raw is ", Cs.cutoff.rw, "..."), quote=FALSE)
print(paste0("Cs cutoff norm is ", Cs.cutoff.norm, "..."), quote=FALSE)

for(chr in chr.v){
  
  load(file=paste0(rawbmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
  rw <- BIN.MX; rm(BIN.MX)
  load(file=paste0(normbmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
  #BIN.MX[BIN.MX==-1] <- 2
  #save(BIN.MX, file=paste0(out.dir, "/", chr, "_", gcb, "_bindata.RData"))
  nrm <- BIN.MX; rm(BIN.MX)
  gc()
  
  if( identical(rw,nrm) ){
    print(paste0(chr, " :BIN.MX weird."), quote=FALSE)
  }

  rw[rw==Cs.cutoff.rw] <- 1
  nrm[nrm==Cs.cutoff.norm] <- 1
  
  if( !identical(rw,nrm) ){
    print(paste0(chr, " :BIN.MX not consistent."), quote=FALSE)
    next
  }
  
  #rm(rw,nrm); gc()
  
  print(paste0(chr, " done!"), quote=FALSE)
  
}

# rm(list=ls()); gc()


