################################################################################
# Function to measure complementarity in terms of matching short k-mers or long 
# span alignment
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(wk.dir, "/lib/getComplementarityAlign.R"))
#- edlib.h - header file
#- edlib.cpp - edlib function source 
#- edlib.NW - Rcpp code calling edlib function
# library(Biostrings)
# library(compiler)
# library(Rcpp)
# sourceCpp(file=paste0(wk.dir, "/lib/edlibNW.cpp"))
# source(paste0(lib.TrantoRextr, "/GEN_getGenomicSeq.R"))  

# source(paste0(wk.dir, "/lib/getComplementarityKmer.R"))
# library(Biostrings)
# library(foreach)
# library(doParallel)
# library(itertools)
# library(compiler)
### FUNCTION ###################################################################
getComplementarityScore <- function(
  contact.mx = contact.mx,
  binlength = 4e4L,
  type = "align", # "kmer" | "align"
  # For type=align, nCPU based on number of chunks
  # For type=kmer, nCPU based on number of contacts 
  nCPU = 1L,
  
  # Alignment-specific arguments
  PATH.genome = "/Users/ltamon/Database/human_genome_unmasked_37.73",
  genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
  fastafile.ending = ".fa",
  silent = FALSE,
  # Number of chunks to split the contacts for alignment
  numChunks = 1,
 
  # K-mer-method-specific arguments
  out.dir = out.dir,
  out.id = paste0(chr, "_Hyb", k, "_", gcb, out.id),
  binkmer.mx = BINKMER.MX,
  # gfree.par file
  gfreeparfile = paste0(gfreepar.dir, "/Gfree_", k, "mer.par"),
  # Kmer length
  kmerlength = 7L,
  saveOut = TRUE
){
  
  if(type=="align"){
    
    print(paste0("Measuring complementarity by alignment in ", numChunks, " chunk/s."),
          quote=FALSE)
    
    SCORE <- getComplementarityAlign(
      PATH.genome=PATH.genome,
      genome.prefix=genome.prefix,
      fastafile.ending=fastafile.ending,
      silent=silent,
      binlength=binlength,
      contact.mx=contact.mx,
      # Number of chunks to split the contacts for alignment
      numChunks=numChunks,
      nCPU=nCPU
    )
    
  } else if(type=="kmer"){
    
    print("Measuring complementarity using kmers...", quote=FALSE)
    
    SCORE <- getComplementarityKmer(
      out.dir=out.dir,
      out.id=out.id,
      contact.mx=contact.mx,
      # Contacts per chromosome
      nCPU=nCPU, 
      binkmer.mx=binkmer.mx,
      # gfree.par file
      gfreeparfile=gfreeparfile,
      # Kmer length
      kmerlength=kmerlength,
      # Length of DNA/Hi-C resolution for Gfree scaling to per-nt
      binlength=binlength,
      saveOut=saveOut
    )
    
    SCORE <- SCORE["NegSumAbsDiff",]
    
  } else {
    stop("Invalid input for type.")
  }
  
  if( length(SCORE)!=nrow(contact.mx) ){
    stop("Checkpoint!")
  } else {
    return(SCORE)
  }
  
}
################################################################################
getComplementarityScore <- cmpfun(getComplementarityScore, options=list(suppressUndefined=TRUE))
################################################################################







