################################################################################
# R API to edlib, lightweight, super fast C/C++ (& Python) library for sequence 
# alignment using edit (Levenshtein) distance
# https://github.com/Martinsos/edlib
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#- edlib.h - header file
#- edlib.cpp - edlib function source 
#- edlib.NW - Rcpp code calling edlib function
# library(Rcpp)
# sourceCpp(file=paste0(wk.dir, "/lib/edlibNW.cpp"))

# library(Biostrings)
# library(foreach)
# library(doParallel)
# library(itertools)
# library(compiler)
# source(paste0(lib.TrantoRextr, "/GEN_getGenomicSeq.R"))  
### FUNCTION ###################################################################
getComplementarityMFE <- function(
  PATH.genome = "/Users/ltamon/Database/human_genome_unmasked_37.73",
  genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
  fastafile.ending = ".fa",
  silent = FALSE,
  binlength = 4e4L,
  contact.mx = contact.mx,
  # Number of chunks to split the contacts for alignment
  numChunks = 15,
  nCPU = 1L,
  RNAduplex.CALL="/home/ltamon/prog/ViennaRNA-2.4.11/bin/RNAduplex",
  filepath.par="./dna_mathews.par",
  id.temp = "chr21"
){
  
  toExport <- c("contact.mx", "binlength")
  if(nCPU > 1){
    registerDoParallel(cores=nCPU)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  # Number of contacts in contact.mx
  ij.len <- nrow(contact.mx)
  
  #### PARALLEL EXECUTION #########
  
  score <- foreach(itr=isplitVector(1:ij.len, chunks=numChunks), .inorder=TRUE,
                   .combine="c", .export=toExport, .noexport=ls()[!ls()%in%toExport]
                   
  ) %op% {
    
    #---------------------------------------
    # Subset of contacts to align
    pair.mx <- contact.mx[itr,]
    pair.len <- nrow(pair.mx)
    
    # Retrieve sequences of the unique bins forming the subset of contacts
    ubins <- unique(c( unique(pair.mx[,1]), unique(pair.mx[,2]) ))
    
    bin.end <- ubins*binlength
    bin.start <- bin.end-binlength+1
    bin.len <- length(bin.end)
    
    ubinsSeq.v <- sapply(X=1:bin.len, simplify=TRUE, FUN=function(i){
      getGenomicSeq(PATH.genome=PATH.genome,
                    genome.prefix=genome.prefix,
                    fastafile.ending=fastafile.ending,
                    chr.id=strsplit(chr,"chr")[[1]][2],
                    remove.other.loads=TRUE,
                    silent=silent, split = FALSE,
                    borders=c(bin.start[i], bin.end[i]) )
    })
    
    print("Sequences obtained.", quote=FALSE)
    
    # Remove chr
    genome.filename <- paste0(genome.prefix, strsplit(chr,"chr")[[1]][2], fastafile.ending)
    eval(parse(text=paste0(
      "rm(", genome.filename, ", bin.start, bin.end);gc()"
    )))
    names(ubinsSeq.v) <- ubins
    
    #---------------------------------------
    # Identify good pairs (no N in contacting sequence)
    goodbins.TF <- !grepl(x=ubinsSeq.v, pattern="N", fixed=TRUE)
    ubinsSeq.v[!goodbins.TF] <- NA
    gpair.ind <- which( pair.mx[,1]%in%ubins[goodbins.TF] & pair.mx[,2]%in%ubins[goodbins.TF] )
    
    rm(genome.filename, ubins, goodbins.TF); gc()
    
    #---------------------------------------
    # Align only the good pairs
    seqi.v <- ubinsSeq.v[ as.character(pair.mx[gpair.ind,1]) ]
    seqj.v <- ubinsSeq.v[ as.character(pair.mx[gpair.ind,2]) ]
    
    # Combine i and j sequences, paired sequences should be next to each other
    seq <- as.vector( rbind(seqi.v, seqj.v) )
    length.v <- rep(binlength, length(seqi.v))
    
    rm(pair.mx, ubinsSeq.v); gc()
    
    #---------------------------------------
    # Main energy calculation (ViennaRNA-RNAduplex) execution
    
    # Since contact information is not strand-specific, calculate
    # two alignment scores: 
    # (Case 1) Between sense strands of seq.i and seq.j
    # (Case 2) Between sense of seq.i and anti-sense of seq.j (equivalent 
    # to alignment of anti-sense of seq.i and sense of seq.j)
    
    #- Case 1: both sense strands
    score.sense <- RNAduplexR(seq=as.vector( rbind(seqi.v, seqj.v) ),
                              RNAduplex.CALL=RNAduplex.CALL,
                              # ID of temp file to avoid overwriting when in parallel
                              id.temp=paste0(id.temp, "_", itr[1]),
                              keep.tmp=FALSE,
                              NuclAc.type="DNA",
                              filepath.par=filepath.par)
                               
    # Reverse complements of j; BOTTLENECK
    seqj.rev.v <- as.vector( reverseComplement(DNAStringSet(seqj.v)) )
    rm(seqj.v)
    
    #- Case 2: sense i - antisense j
    score.rev <- RNAduplexR(seq=as.vector( rbind(seqi.v, seqj.rev.v) ),
                            RNAduplex.CALL=RNAduplex.CALL,
                            id.temp=paste0(id.temp, "_", itr[1]),
                            keep.tmp=FALSE,
                            NuclAc.type="DNA",
                            filepath.par=filepath.par)
                              
    rm(seqi.v, seqj.rev.v); gc()
    
    # Return the more negative value (minimum energy) between two cases; 
    # bad pairs denoted by NA
    score.chunk <- rep(NA, pair.len)
    score.chunk[gpair.ind] <- apply(X=cbind(score.sense, score.rev), MARGIN=1, 
                                    FUN=min, na.rm=TRUE)
    
    print("Chunk done!", quote=FALSE)
    
    ## Normalise scores to bin length
    #return(score.chunk/binlength)
    return(score.chunk)
    
    #---------------------------------------
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  return(score)

}
################################################################################
getComplementarityMFE <- cmpfun(getComplementarityAlign, options=list(suppressUndefined=TRUE))
################################################################################

# rm(list=ls())
 