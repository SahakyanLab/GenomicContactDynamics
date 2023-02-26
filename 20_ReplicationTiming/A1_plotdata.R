################################################################################
# Generate vals object with rt calculations per chr. Lengths of vals object for a 
# chromosome should all be equal to the total chr bins based on HiC resolution. 
# The code also expects that the resolution of RT data should be equal to HiC 
# resolution. Code checks that there are no rt bins exceeding the number of chr
# bins based on resolution.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/20_ReplicationTiming")
out.dir = paste0(wk.dir, "/out_plotdata")

filter.id = "mincelllineWithData59"
reptimePath = paste0(wk.dir, "/out_transform_rtdata/RT_data_hg19_", filter.id, ".RData")
chrLen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
rt.calctype = c("mean.all", "mean.tumor", "mean.nontumor",
                "median.all", "median.tumor", "median.nontumor")
chrs = paste0("chr", c(1:22, "X"))
nCPU = 5 # Number of combinations between chrs and rt.calctype
HiC.res = 40000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
#source(paste0(wk.dir, "/lib/checkANDfilterRTdata.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
rt <- readRDS(file=reptimePath)
rt$chroms <- as.character(rt$chroms)
rt$bin <- ceiling( (rt$starts + HiC.res - 1) / HiC.res )

## Check and filter RT data
#rt <- checkANDfilterRTdata(rt=rt, HiC.res=HiC.res, chrs=chrs, chrLen.file=chrLen.file)

# Generate .RData per chr and per rt calc

combis.df <- expand.grid(chr=chrs, rt.calctype=rt.calctype, stringsAsFactors=F)
combis.len <- length(combis.df[,1])

chrLen.df <- read.delim(file=chrLen.file, header=T)  
out.id <- paste0(filter.id, "_", HiC.res, "bpHiCres")

toExport <- c("combis.df", "rt", "out.dir", "out.id", "chrLen.df") 
#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:combis.len, chunks=nCPU), .inorder=F,
        .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in itr){
    
    rt.calctype <- combis.df$rt.calctype[[i]]
    rt.calc <- strsplit(rt.calctype, split=".", fixed=T)[[1]][1]
    rt.type <- strsplit(rt.calctype, split=".", fixed=T)[[1]][2]
    
    #
    
    chr <- combis.df$chr[[i]]
    chr.tot.bin <- chrLen.df$bins.40kb[chrLen.df$chromosome == chr]
    is.rt.chr <- rt$chroms == chr
    
    #
    
    if( max(rt$bin[is.rt.chr]) > chr.tot.bin ){
      
      rm(rt)
      stop(paste0(chr, " ", rt.calctype, ": chr rt max bin > chr total bins."))
      
    }
    
    #
    
    vals <- rep(NA, times=chr.tot.bin)
    vals[ rt$bin[is.rt.chr] ] <- rt[[rt.calctype]][is.rt.chr]
    
    save(vals, file=paste0(out.dir, "/", chr, "_", rt.type, "_", rt.calc, "_", out.id, ".RData"))
    
    rm(i, chr, rt.calc, vals)
    
  }
  
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()