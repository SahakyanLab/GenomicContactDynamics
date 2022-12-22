################################################################################
# Get Cp data overlapping with each converted hg38 to hg19 contacts from 
# F1_liftover.R
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) 
options(warn=1)

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
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity")
ijLO.dir = paste0(wk.dir, "/out_liftover") # 1-based
out.dir = paste0(wk.dir, "/out_assignCp")
### OTHER SETTINGS #############################################################
src.id = "hg38ToHg19_LOwidth.min.bp30000_kmer_min2Mb"
chrs = "chrCHRREPLACE"
coord.system = 1 
# PERSIST.MX
gcb = "min2Mb"
bin.len = 40000 

# Arguments dependent on cons.methd
cons.methd = 2 # i.e. data.table method of getConsensusCpOfPairs()
consensus.FUN = function(z){paste(z, collapse=";")} 
consensus.FUN.id = paste0("cons.methd", cons.methd, "_finalStringOfCps") 
nCPU = 5 # Contacts
numChunks = 500 # Process p.coord.mx in chunks 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/getConsensusCpOfPairs.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chrs){
  
  # Load PERSIST.MX
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  
  ij.mx <- cbind(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
  rm(PERSIST.MX)
  
  # Load p.coord.mx
  load(paste0(ijLO.dir, "/", chr, "_", src.id, ".RData"))
  ijLO.mx <- data.matrix(ijLO.df[,-1])
  #ijLO.mx <- ijLO.mx[1:1000,] # REMOVE
  rm(ijLO.df)

  # Get Cps overlapping with each region pairs in p.coord.mx
  
  #### PARALLEL EXECUTION #########
  
  toExport <- c("bin.len", "ij.mx", "coord.system", "ijLO.mx", "consensus.FUN", "cons.methd", "nCPU")
  p.coord.len <- length(ijLO.mx[,1]) 
  
  consensusCp <- foreach(itr=isplitVector(1:p.coord.len, chunks=numChunks), .inorder=T,
                         .combine="rbind", .export=toExport, .noexport=ls()[!ls()%in%toExport]
                            
  ) %op% {
    
    print(paste0("Chunk starting index ", itr[[1]], "..."), quote=F)
    chunk <- getConsensusCpOfPairs(bin.len=bin.len,
                                   ij.mx=ij.mx, # [i, j, Cp] 
                                   coord.system=coord.system,
                                   p.coord.mx=ijLO.mx[itr,], # REMOVE # 1-based CS, [start.a, end.a, start.b, end.b]
                                   consensus.FUN=consensus.FUN, # Depends on method, see description above
                                   methd=cons.methd, # INT, 1 or 2
                                   nCPU=nCPU)
    chunk$pair.ind <- itr
    
    return(chunk)
                                    
  }
  
  ### END OF PARALLEL EXECUTION ###
    
  ij.pair.len <- length(ijLO.mx[,1]) # REMOVE
  if( !identical(as.integer(consensusCp[,"pair.ind"]), 1:ij.pair.len) ){
    rm(consensusCp)
    stop(paste(chr, ": consensusCp order not consistent with ijLO.df/mx."))
  }

  consensusCp <- consensusCp[,setdiff(colnames(consensusCp), "pair.ind"), drop=F]
  save(consensusCp, file=paste0(out.dir, "/", chr, "_", src.id, "_", consensus.FUN.id, "Cp.RData"))
  
  print(paste0(chr, ": done!"), quote=F)
  
  rm(ijLO.mx, consensusCp, chr)
  gc()
  
}

# rm(list=ls()); gc()
