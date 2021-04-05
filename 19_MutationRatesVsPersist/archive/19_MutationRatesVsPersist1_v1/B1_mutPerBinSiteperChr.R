################################################################################
# Mutation frequency per nucleotide site in a bin (done per chr)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
src.dir = paste0(wk.dir, "/out_filter")
out.dir = paste0(wk.dir, "/out_mutPerBinSiteperChr")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
src.id = "hg38ToHg19"
gcb = "min2Mb"
chr.v = paste("chr", 1:22, sep="")
bin.len = 40000
nCPU = 4L # Number of bins
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
library(IRanges)
library(GenomicRanges)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(wk.dir, "/lib/olapPerGenomicPos.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrLen.df <- read.delim(file=chrLenfile, header=TRUE)

#-------------------Collect contact bins from all chr
BIN.MX <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
  # Load PERSIST.MX
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  ubins <- unique( c(unique(PERSIST.MX$hits$i), unique(PERSIST.MX$hits$j)) )
  rm(PERSIST.MX)
  data.frame(chr=chr, bin=ubins, stringsAsFactors=FALSE, row.names=NULL)
})
BIN.MX <- do.call("rbind", BIN.MX)
rownames(BIN.MX) <- NULL

#-------------------Feature/mutation file
load(file=paste0(src.dir, "/CosmicNCV_", src.id, "_final_2000.RData"))
#load(file=paste0(src.dir, "/CosmicNCV_", src.id, "_final.RData"))
ncv.df <- ncv.df[,c("ID_SAMPLE", "chr", "start", "end", "MUT")]

# Filter ncv.df by mutation type
mut.v <- c("All", unique(ncv.df$MUT))
if( is.null(chr.v) ){ chr.v <- unique(ncv.df$chr) }
ncv.df <- ncv.df[,c("chr", "start", "end", "MUT")]

#-------------------Measure mutation frequency per base in each 40-kb bin
bin.half <- bin.len/2
pos.v <- (-bin.half+1L):bin.half

for(mut in mut.v){
  
  if(mut=="All"){
    mut.TF <- rep(TRUE, times=nrow(ncv.df))
  } else if( mut%in%c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T") ){
    mut.TF <- ncv.df$MUT==mut
  } else {
    stop("Invalid mutation type notation.")
  }
  
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=TRUE)
  
  for(chr in chr.v){
    
    maxPos <- as.integer(chrLen.df[chrLen.df[,"chromosome"]==chr, "length.bp"])
    
    ubins <- BIN.MX[BIN.MX$chr==chr,"bin"]
    bin.end <- bin.start <- (ubins*bin.len)-(bin.half)
  
    incl.TF <- mut.TF & ncv.df$chr==chr
    
    if( sum(incl.TF)==0L ){
      
      print(paste0(chr, " skipped"), quote=FALSE)
      next
      
    } else {
      
      FETA.MX <- olapPerGenomicPos(
        # Feature of interest
        foi.id=paste0(chr, "_", mut.id),
        foi.df=ncv.df[incl.TF,c("chr", "start", "end")],
        start=bin.start,
        end=bin.end,
        res=1,
        pos.v=pos.v,
        nCPU=nCPU,
        minPos=1L,
        maxPos=maxPos,
        # Overlap parameters
        min.olap=1L,
        # This default is fixed regardless of os for R.3.5.0 and up
        max.gap=-1L, 
        type.olap="any"
      )
      
    }
    
    rm(incl.TF, bin.start, bin.end)
    
    FETA.MX <- cbind(bin=ubins, FETA.MX)
    dimnames(FETA.MX)[[2]] <- c("bin", as.character(format(x=pos.v, scipen=FALSE)))
    
    save(FETA.MX, file=paste0(out.dir, "/temp/", chr, "_", gcb, "_", mut.id, ".RData"))
    
    rm(FETA.MX); gc()
    
  } # chr.v for loop end
  
  print(paste0(mut.id, " done!"), quote=FALSE)
  rm(mut.TF, mut.id)
  
} # mut.v for loop end

# rm(list=ls()); gc()


