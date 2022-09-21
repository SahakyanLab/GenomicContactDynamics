################################################################################
# Per chr, return complementarity measures per contact
### FUNCTION ###################################################################
feat_complement <- function(gcb, chr, PERSIST.MX, CII.dir){
  
  contact.mx <- data.matrix(PERSIST.MX$hits[,c("i", "j")])
  rm(PERSIST.MX)
  
  #
  load(paste0(CII.dir, "/", chr, "_kmer_", gcb, ".RData"))
  
  CII.ind <- match(rownames(contact.mx), table=rownames(CII.MX))
  if( !identical(as.numeric(contact.mx), 
                 as.numeric(CII.MX[CII.ind,c("i", "j")])) ){
    rm(CII.MX); stop("feat_complement(): Checkpoint CIIkmer.")
  }
  
  feat.mx <- CII.MX[CII.ind, c("C||", "Gfree")]
  dimnames(feat.mx)[[2]] <- c("CIIkmer", "CIIG")
  rm(CII.MX)
  
  #
  load(paste0(CII.dir, "/", chr, "_align_", gcb, ".RData"))
  
  if( !identical(as.numeric(contact.mx), 
                 as.numeric(CII.MX[CII.ind,c("i", "j")])) ){
    rm(CII.MX); stop("feat_complement(): Checkpoint CIIalign.")
  }
  
  feat.mx <- cbind(feat.mx, CIIalign=CII.MX[CII.ind,"C||"])
  dimnames(feat.mx)[[2]] <- paste0("ij_", dimnames(feat.mx)[[2]])
  
  return(feat.mx)

}

# rm(list=ls()); gc()
