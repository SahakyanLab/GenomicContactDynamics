################################################################################
# Combine CII.MX (type=kmer) and HYB.MX; normalise align-based scores to length
# of region; make 0.5Mb version of CII.MX by modifying the Cp column
# deva, R/3.6.0-newgcc, gcc/4.9.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    #wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    data.dir =  "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    #wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    data.dir =  "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
compl.dir = paste0(wk.dir, "/out_constraints")
hyb.dir = paste0(wk.dir, "/out_constraints")
out.dir = paste0(wk.dir, "/out_constraints/merged_final")
#persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
### OTHER SETTINGS #############################################################
chr.v = "chr17" #paste("chr", c(1:16, 18:22, "X"), sep="")
exclude = c(`chr17`=2030)
gcb = "min2Mb"
makeothergcb = TRUE
type = "align" # kmer | align
kmer.len = 7L
bin.len = 40000
affix = "" #"_ijShuffled" 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
gcb.other <- c("min2Mb", "min05Mb")[c("min2Mb", "min05Mb")!=gcb]
for(chr in chr.v){
  # CII.MX
  load(paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"))
  if(type=="kmer"){
    # HYB.MX
    load(paste0(hyb.dir, "/", chr, "_Hyb", kmer.len, "_", gcb, affix, ".RData"))
    if(nrow(CII.MX)!=ncol(HYB.MX)){ stop("Length of CII.MX and HYB.MX different.") }
    if( !identical(as.numeric(CII.MX[,4]),as.numeric(HYB.MX[3,])) ){
      stop("Complementarity values of CII.MX and HYB.MX different.")
    }
    CII.MX <- cbind(CII.MX, t(HYB.MX[rownames(HYB.MX)!="NegSumAbsDiff",])
    )
    rm(HYB.MX); gc()
    
  } else if(type=="align"){
    # Normalise to length of region and exclude contacts with different lengths
    # of regions
    if( chr%in%names(exclude) ){
      exclij.TF <- CII.MX[,"i"]==exclude[chr] | CII.MX[,"j"]==exclude[chr]
      CII.MX[exclij.TF,"C||"] <- NA
      print(paste0(chr, ": ", sum(exclij.TF),  " contact/s formed by bin ", 
                   exclude[chr], " excluded."))
    }
    CII.MX[,"C||"] <- CII.MX[,"C||"]/bin.len
  } else {
    stop("Type can only be 'kmer' or 'align'")
  }
  save(CII.MX, file=paste0(out.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"))
  #---------------------------------------
  if(makeothergcb){
    print(paste0("Making ", gcb.other, " version..."), quote=FALSE)
    # PERSIST.MX min05Mb to add Cp min05
    load(paste0(persist.dir, "/", chr, "_Persist_", gcb.other, affix, ".RData"))
    Cp.v <- PERSIST.MX$ntis
    names(Cp.v) <- rownames(PERSIST.MX$hits)
    rm(PERSIST.MX); gc()
    CII.MX[names(Cp.v),"Cp"] <- Cp.v
    if( length(Cp.v)!=sum(!is.na(CII.MX[,3])) ){
      stop(paste0(gcb.other, ": Number of contacts with cp different."))
    } 
    save(CII.MX, file=paste0(out.dir, "/", chr, "_", type, "_", gcb.other, affix, ".RData"))
  }
  
  print(paste0(chr, " done!"), quote=FALSE)
}

# rm(list=ls())






