################################################################################
# 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" 
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/27_EvolutionaryAnalyses")
CII.state.dir = paste0(wk.dir, "/out_complementarityVsContactValue") 
CII.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity/out_constraints_hg38_GfreeSingleNorm/merged_final")
cp.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity/out_consensusCp")
cf.dir = paste0(data.dir, "/Phylo-HMRF/out_extractFromHic_alldata")
out.dir = paste0(wk.dir, "/out_complementarityVsContactValue_unabridged") 
### OTHER SETTINGS #############################################################
chrs = paste0("chr", c("X"))
CII.state.id = "min2Mb_genome_state_Phylo-HMRF_contact50K_norm_KR"
# If no CII state data (VAL.MX), get complementarity values from CII.MX
gcb = "min2Mb"
compl.data.ids = "kmer"
compl.types = c("C||", "Gfree", "sdDifference")
cp.id = "hg38ToHg19_LOwidth.min.bp30000_kmer_min2Mb_consensusCp"
cf.id = "contact50K_NONE_KR"
out.id = paste0(CII.state.id, "_CP", cp.id, "_CF", cf.id)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cf.species = c("hg38.downsample", "GSM3685696_gorGor4_mapping",
               "GSM3685695_panPan2_mapping", "GSM3685694_panTro5_mapping")
cf.nmes = c("norm", "KR")
tmp <- expand.grid(cf.species, cf.nmes)
cf.cols.append <- paste0(tmp$Var1, "_", tmp$Var2)

for(chr in chrs){
  
  CII.state.fle <- paste0(CII.state.dir, "/", chr, "_", CII.state.id, ".RData")
  if( file.exists(CII.state.fle) ){
    load(CII.state.fle)
  } else {
    
    VAL.MX <- list()
    for(id in compl.data.ids){
      
      load(paste0(CII.dir, "/", chr, "_", id, "_", gcb, ".RData"))
      VAL.MX[["i"]] <- CII.MX[,"i", drop=F]
      VAL.MX[["j"]] <- CII.MX[,"j", drop=F]
      dimnames(CII.MX)[[2]] <- paste0(id, ".", dimnames(CII.MX)[[2]])
      VAL.MX[[id]] <- CII.MX[ ,intersect(paste0(id, ".", compl.types), dimnames(CII.MX)[[2]]), drop=F ] 
      rm(CII.MX)
      
      print(paste0(chr, " ", id, ": Adding data."), quote=F)
      
    }
    names(VAL.MX) <- NULL
    VAL.MX <- do.call("cbind.data.frame", VAL.MX)
    VAL.MX.append <- matrix(data=NA, nrow=length(VAL.MX[,1]), ncol=length(cf.species) * 2 + 1,
                            dimnames=list(NULL, c("PHYLOHMRFstate", cf.cols.append)))
    VAL.MX <- cbind(VAL.MX, VAL.MX.append)
    rm(VAL.MX.append)
  
  }
  
  ijlen.valmx <- length(VAL.MX[,1])
  collen.valmx <- length(VAL.MX[1,])
  
  for(cf.spec in cf.species){
    
    load(paste0(cf.dir, "/", chr, "_", cf.spec, "_", cf.id, ".RData"))
    
    if( !identical(ijlen.valmx, length(HIC.ALLIJ[,1])) ){
      rm(HIC.ALLIJ)
      stop(paste0(chr, ": Different HIC.ALLIJ and VAL.MX row count i.e. 
                  count of contacts."))
    }
    
    cf.spec.id <- paste0(cf.spec, "_KR")
    VAL.MX[[cf.spec.id]] <- HIC.ALLIJ[[cf.spec.id]]
    
    cf.spec.id <- paste0(cf.spec, "_NONE")
    VAL.MX[[cf.spec.id]] <- HIC.ALLIJ[[cf.spec.id]]
    
    rm(HIC.ALLIJ, cf.spec.id)
    
    print(paste0(chr, " ", cf.spec, ": HIC.ALLIJ data added."), quote=F)
    
  }
  
  #if( !identical(collen.valmx, length(VAL.MX[1,])) ){
  #  rm(VAL.MX)
  #  stop(paste0(chr, ": Addition of KR data made new column/s."))
  #}
  
  # Add Cp data
  
  load(paste0(cp.dir, "/", chr, "_", cp.id, ".RData"))
  
  ## Check if new consensus Cp data matches with old consensus Cp data
  #is_withstate <- !is.na(VAL.MX$PHYLOHMRFstate)
  #old.cp.nmes <- colnames(VAL.MX)[grepl(colnames(VAL.MX), pattern=".consCp", fixed=T)]
  #check.TF <- sapply(X=old.cp.nmes, simplify=T, FUN=function(nme){
  #  !identical(VAL.MX[[nme]][is_withstate], consensusCp[,nme][is_withstate])
  #})
  
  #if( any(unname(check.TF)) ){
  #  rm(VAL.MX)
  #  stop(paste0(chr, ": Old and new (final) consensus Cp data not matching."))
  #}
  
  #old.cp.nmes <- colnames(VAL.MX)[grepl(colnames(VAL.MX), pattern=".consCp", fixed=T)]
  #VAL.MX <- VAL.MX[, setdiff(colnames(VAL.MX), old.cp.nmes) ]  
  
  VAL.MX <- cbind(VAL.MX, consensusCp)
  rm(consensusCp)
  #rm(consensusCp, is_withstate, check.TF, old.cp.nmes)
  
  print(paste0(chr, ": Final consensus Cp data added."), quote=F)
  
  #
  
  save(VAL.MX, file=paste0(out.dir, "/", chr, "_", out.id, ".RData"))
  rm(VAL.MX)
  
  print(paste0(chr, " done!"), quote=F)
  
}

# rm(list=ls()); gc()