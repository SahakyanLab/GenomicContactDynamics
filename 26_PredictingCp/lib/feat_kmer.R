################################################################################
# Per chr, return kmer counts per contact
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(foreach)
# library(doParallel)
# library(itertools)
# library(Biostrings)
# library(compiler)
# source(paste0(lib, "/makeKmerStrandInvar.R"))
### FUNCTION ###################################################################
feat_kmer <- function(BINKMER.MX,
                      contact.mx = 'i-j 2-column matrix',
                      kmer.len = 'kmer length'
                      ){
  
  BINKMER.MX <- makeKmerStrandInvar(mx=BINKMER.MX, nCPU=1, removeRevcomp=T)
  
  kmer.vec <- setdiff( dimnames(BINKMER.MX)[[2]], c("bins", "startpos", "endpos", "numUMChar") )
  
  # Checks 
  kmer.num.noRC <- (4^kmer.len) / 2
  is_wrong <- c(any(duplicated(dimnames(BINKMER.MX)[[2]])), 
                length(kmer.vec) != kmer.num.noRC,
                unique(nchar(kmer.vec)) != kmer.len,
                any(!is.numeric(contact.mx[,1])),
                any(!is.numeric(contact.mx[,2])),
                any(!is.numeric(BINKMER.MX[,"bins"]))
                )
  if( any(is_wrong) ){
    rm(BINKMER.MX)
    stop("feat_kmer(): Unexpected behaviours based on kmer length.")
  }
  
  #
  i.binkmerInd <- match( contact.mx[,1], table=BINKMER.MX[,"bins"] )
  j.binkmerInd <- match( contact.mx[,2], table=BINKMER.MX[,"bins"] )
  
  BINKMER.MX <- BINKMER.MX[,kmer.vec]
  feat.mx <- cbind(BINKMER.MX[i.binkmerInd, ], BINKMER.MX[j.binkmerInd, ])
  dimnames(feat.mx) <- list(dimnames(contact.mx)[[1]], 
                            c(paste0("i_", kmer.vec), paste0("j_", kmer.vec)))
  
  return(feat.mx)

}

# rm(list=ls()); gc()
