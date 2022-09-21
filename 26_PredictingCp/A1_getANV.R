################################################################################
# Get accumulated natural vector (ANV) (doi: 10.3389/fgene.2019.00234) for all 
# 40-kb regions of human genome.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/26_PredictingCp")
genome.dir = paste0(data.dir, "/human_genome_unmasked_37.73")
out.dir = paste0(wk.dir, "/out_getANV")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(22:1, "X"))
nCPU = 1
bin.len = 40000
genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome."
fastafile.ending = ".fa"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
library(Biostrings)
source(paste0(lib, "/TrantoRextr/UTIL_readLinesFast.R"))   
source(paste0(lib, "/TrantoRextr/GEN_readfasta.R"))  
source(paste0(lib, "/TrantoRextr/GEN_loadGenome.R"))                             
source(paste0(lib, "/TrantoRextr/GEN_getGenomicSeq.R"))  
library(Rcpp)
sourceCpp(file=paste0(wk.dir, "/lib/anv_encoding.cpp"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrlen.df <- read.delim(chrlen.file, header=T, stringsAsFactors=F)
chrlen.df$tot.bin <- ceiling(chrlen.df$length.bp / bin.len)

if( !identical(as.numeric(chrlen.df$bins.40kb), as.numeric(chrlen.df$tot.bin)) ){
  stop("Checkpoint 1.")
  rm(chr.v)
}
chr.v.len <- length(chr.v)

bases <- c("A", "C", "G", "T")
combi <- unlist(lapply(X=combn(x=bases, m=2, simplify=F), paste, collapse=""))

toExport <- c("chr.v", "chrlen.df", "bin.len", "genome.dir", "genome.prefix", 
              "fastafile.ending", "bases", "combi")

#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), .inorder=F, 
        .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {

  for(i in itr){
    
    chr <- chr.v[i]
    tot.bin <- chrlen.df$tot.bin[chrlen.df$chromosome == chr]
    ubins <- 1:tot.bin
    bin.end <- ubins * bin.len
    bin.start <- bin.end - bin.len + 1
    bin.end[tot.bin] <- chrlen.df$length.bp[chrlen.df$chromosome == chr]
    
    ubinSeqs <- sapply(X=ubins, simplify=T, FUN=function(ind){
      getGenomicSeq(PATH.genome=genome.dir,
                    genome.prefix=genome.prefix,
                    fastafile.ending=fastafile.ending,
                    chr.id=strsplit(chr,"chr")[[1]][2],
                    remove.other.loads=T, silent=T, split=F,
                    borders=c(bin.start[ind], bin.end[ind]) )
    })
    
    is_validSeq <- !grepl(ubinSeqs, pattern="N")
    
    anv.mx <- mapply(FUN=rcpp_anv_encode, `char_string`=ubinSeqs[is_validSeq], SIMPLIFY=T, USE.NAMES=F)
    
    rm(ubinSeqs)
    
    ANV.MX <- matrix(NA, nrow=tot.bin, ncol=18)
    ANV.MX[is_validSeq,] <- t(anv.mx)
    
    rm(anv.mx, is_validSeq)
  
    dimnames(ANV.MX) <- list( as.character(ubins),
                              c(paste0("n", bases), paste0("zeta", bases), paste0("D", bases), paste0("cov", combi)) )
    ANV.MX <- cbind(bins=ubins, startpos=bin.start, endpos=bin.end, numUMChar=(bin.end - bin.start + 1), ANV.MX)
    
    save(ANV.MX, file=paste0(out.dir, "/", chr, "_ANV.RData"))
    
    rm(ANV.MX, ubins, bin.start, bin.end, chr, tot.bin)
    
  }
  
}
### END OF PARALLEL EXECUTION ###

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()
