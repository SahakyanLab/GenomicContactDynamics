################################################################################
# Per chr, return table for ML
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(foreach)
# library(doParallel)
# library(itertools)
# library(Biostrings)
# library(compiler)
# source(paste0(lib, "/makeKmerStrandInvar.R"))
# source(paste0(wk.dir, "/lib/feat_kmer.R"))
# source(paste0(wk.dir, "/lib/feat_complement.R"))
# source(paste0(wk.dir, "/lib/feat_anv.R"))
### FUNCTION ###################################################################
feat_makeTablePerChr <- function(anv.dir, binkmer3.dir, binkmer1.dir, CII.dir, 
                                 persist.dir, chr, gcb
                                 ){
  
  TBL <- list()
  
  load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  
  TBL$LABEL <- PERSIST.MX$ntis
    
  #
  TBL$grp.compl <- feat_complement(gcb=gcb, chr=chr, PERSIST.MX=PERSIST.MX, CII.dir=CII.dir)
  
  #
  contact.mx <- data.matrix(PERSIST.MX$hits[,c("i", "j")])
  rm(PERSIST.MX)
  
  #
  euclidean <- function(a, b){ sqrt(sum((a - b)^2)) }
  load(paste0(anv.dir, "/", chr, "_ANV.RData"))
  TBL$grp.anv <- feat_anv(ANV.MX=ANV.MX, contact.mx=contact.mx, distanceOnly=F, distFUN=euclidean)
  
  TBL$grp.compl <- cbind(ij_anvdist=TBL$grp.anv[,"ij_anvdist"], TBL$grp.compl)
  TBL$grp.anv <- TBL$grp.anv[,-(dimnames(TBL$grp.anv)[[2]] == "ij_anvdist")]
           
  #
  load(paste0(binkmer3.dir, "/", chr, "_BinKmer3.RData"))
  TBL$grp.kmer3 <- feat_kmer(BINKMER.MX=BINKMER.MX, contact.mx=contact.mx, kmer.len=3)
  
  #
  load(paste0(binkmer1.dir, "/", chr, "_BinKmer1.RData"))
  TBL$grp.kmer1 <- feat_kmer(BINKMER.MX=BINKMER.MX, contact.mx=contact.mx, kmer.len=1)
  
  # GC content
  TBL$grp.GC <- cbind(i_GC=TBL$grp.kmer1[,intersect(colnames(TBL$grp.kmer1), c("i_G", "i_C"))] / 
                        rowSums(TBL$grp.kmer1[,grepl("i_", x=colnames(TBL$grp.kmer1), fixed=T)], na.rm=F),
                      j_GC=TBL$grp.kmer1[,intersect(colnames(TBL$grp.kmer1), c("j_G", "j_C"))] / 
                        rowSums(TBL$grp.kmer1[,grepl("j_", x=colnames(TBL$grp.kmer1), fixed=T)], na.rm=F))
  TBL$grp.GC <- cbind(TBL$grp.GC, ij_meanGC=rowMeans(TBL$grp.GC, na.rm=F))
  
  TBL <- do.call("cbind.data.frame", list(TBL, stringsAsFactors=F))
  
  #  Checks
  is_wrong <- c(any(duplicated(colnames(TBL))))
  
  if( any(is_wrong) ){
    rm(TBL); stop("feat_makeTablePerChr(): Checkpoint.")
  } else {
    return(TBL)
  }
  
}

# rm(list=ls()); gc()
