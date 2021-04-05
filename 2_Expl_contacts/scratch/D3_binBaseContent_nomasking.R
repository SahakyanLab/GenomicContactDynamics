################################################################################
# Calculate base content per bin
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/2_Expl_contacts"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_Expl_contacts"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
out.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binBaseContent")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(1:22, "X"))
nCPU = 2L
genomePath = paste0(data.dir, "/human_genome_unmasked_37.73")
genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome."
fastafile.ending = ".fa"
bin.len = 4e4
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/TrantoRextr/GEN_loadGenome.R"))                             #
source(paste0(lib, "/TrantoRextr/GEN_readfasta.R"))                              #
source(paste0(lib, "/TrantoRextr/UTIL_readLinesFast.R"))   
source(paste0(lib, "/TrantoRextr/GEN_getGenomicSeq.R"))                          #
source(paste0(lib, "/TrantoRextr/GEN_getKmers.R")) 
source(paste0(lib, "/TrantoRextr/GEN_getKmerCountsPerInterval.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrLen.df <- read.delim(file=chrLenfile, header=TRUE)

chr.v.len <- length(chr.v)

toExport <- c("chr.v", "chrLen.df", "bin.len", "genomePath", "genome.prefix", 
              "fastafile.ending", "out.dir")

#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
        .inorder=FALSE, .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in itr){
    
    chr <- chr.v[i]
    totlen <- chrLen.df$length.bp[chrLen.df$chromosome==chr]
    totbin <- ceiling(totlen/bin.len)
    if( totbin!=chrLen.df$bins.40kb[chrLen.df$chromosome==chr] ){
      stop("Checkpoint 1.")
    }
    
    endpos <- (1:totbin)*bin.len
    startpos <- endpos-bin.len+1L
    endpos[totbin] <- totlen
    lastbin.len <- endpos[totbin]-startpos[totbin] + 1
    
    BINKMER.MX <- getKmerCountsPerInterval(
      chr=chr,
      startpos=startpos,
      endpos=endpos,
      K=1,
      PATH.genome=genomePath,
      genome.prefix=genome.prefix,
      fastafile.ending=fastafile.ending,
      silent=TRUE
    )
    rsum <- rowSums(x=BINKMER.MX[,-totbin], na.rm=FALSE)
    if( !all(rsum%in%c(NA, 40000, lastbin.len)) ){
      stop(paste0(chr, ": Bin length not 40000."))
    }
    
    rsum[totbin] <- lastbin.len
    rsum[is.na(rsum)] <- bin.len
    
    BINKMER.MX <- cbind(bins=1:totbin, startpos=startpos, endpos=endpos,
                        numUMChar=rsum, BINKMER.MX)
    save(BINKMER.MX, file=paste0(out.dir, "/", chr, "_BinKmer1.RData"))
    
    rm(rsum, BINKMER.MX, startpos, endpos, totbin, totlen); gc()
    
  }
  
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()