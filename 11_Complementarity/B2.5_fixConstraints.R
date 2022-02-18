################################################################################
# Combine CII.MX (type=kmer) and HYB.MX; normalise align-based scores to 
# Hi-C resolution (kmer-based ones were already normalised); set values involving 
# last bin (using exclude argument) to NA because the length of the last bin is 
# shorter than the Hi-C resolution so the normalisation will be faulty (originally 
# the code does this only for align-based scores but I realised when doing ath 
# that this also applies for k-mer based ones.); for human datasets, make 0.5Mb 
# version of CII.MX by modifying the Cp column.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/11_Complementarity"
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
compl.dir = paste0(wk.dir, "/z_ignore_git/out_constraints_ath_20kb")
hyb.dir = paste0(wk.dir, "/z_ignore_git/out_constraints_ath_20kb")
out.dir = paste0(wk.dir, "/z_ignore_git/out_constraints_ath_20kb/merged_final")
#persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
persist.dir = paste0(data.dir, "/HiC_features_arabidopsis_HiC_NORMpc/features_20kb")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", 1:5) #paste("chr", c(1:16, 18:22, "X"), sep="")
exclude = setNames(object=c(1522, 985, 1173, 930, 1349),
                   nm=paste0("chr", 1:5))
gcb = "min0Mb"
possible.gcb = "min0Mb"
makeothergcb = FALSE
type = "align" # kmer | align
kmer.len = 7L
bin.len = 20000
affix = "" #"_ijShuffled" 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Applicable for human data
gcb.other <- possible.gcb[possible.gcb!=gcb]

for(chr in chr.v){
  
  # CII.MX
  load(paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"))
  
  if(type=="kmer"){
    
    # HYB.MX
    load(paste0(hyb.dir, "/", chr, "_Hyb", kmer.len, "_", gcb, affix, ".RData"))
    if( nrow(CII.MX)!=ncol(HYB.MX) ){ stop("Length of CII.MX and HYB.MX different.") }
    if( !identical(as.numeric(CII.MX[,4]),as.numeric(HYB.MX[3,])) ){
      stop("Complementarity values of CII.MX and HYB.MX different.")
    }
    CII.MX <- cbind(CII.MX, t(HYB.MX[rownames(HYB.MX)!="NegSumAbsDiff",])
    )
    rm(HYB.MX); gc()
    
  } else if(type=="align"){
    
    CII.MX[,"C||"] <- CII.MX[,"C||"]/bin.len
    print("Divided C||align values by bin.len.")
    
  } else {
    stop("Type can only be 'kmer' or 'align'")
  }
  
  # Exclude contacts with bins differing in length
  if( chr%in%names(exclude) ){
    
    exclij.TF <- CII.MX[,"i"]==exclude[[chr]] | CII.MX[,"j"]==exclude[[chr]]
    CII.MX[exclij.TF,"C||"] <- NA
    print(paste0(chr, ": ", sum(exclij.TF),  " contact/s formed by bin ", 
                 exclude[chr], " excluded."))
    
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

# rm(list=ls()); gc()






