################################################################################
# Generate vals object with rt calculations per chr. Lengths of vals object for a 
# chromosome should all be equal to the total chr bins based on HiC resolution. 
# The code also expects that the resolution of RT data should be equal to HiC 
# resolution. Before generating .vals object, all values of bins/regions 
# not satisfying filtering criteria are set to NA ( see checkANDfilterRTdata() ).
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
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/20_ReplicationTiming")
out.dir = paste0(wk.dir, "/out_plotdata")

rt.type = "all" # "all" | "tumor" | "nontumor"
reptimePath = paste0(data.dir, "/replication_timing/out_clustering_combined/hg19/", 
                     rt.type, "/RT_data_hg19.RData")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
rt.calcs = c("mean", "median", "sd", "set.count", "point.count")
chrs = paste0("chr", c(21:22))
nCPU = 1 # Number of combinations between chrs and rt.calcs
HiC.res = 40000
Cp.v = 1:21
filter.id = "normNot1_setpointcountGrEq3"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(wk.dir, "/lib/checkANDfilterRTdata.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
rt <- readRDS(file=reptimePath)

# Check and filter RT data
rt <- checkANDfilterRTdata(rt=rt, HiC.res=HiC.res, chrs=chrs, chrLenfile=chrLenfile)

# Generate .RData per chr and per rt calc

combis.df <- expand.grid(chr=chrs, rt.calc=rt.calcs, stringsAsFactors=F)
combis.len <- length(combis.df[,1])
  
out.id <- paste0(filter.id, "_", HiC.res, "bpHiCres")
toExport <- c("combis.df", "rt", "out.dir", "out.id")
#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:combis.len, chunks=nCPU), .inorder=F,
        .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in itr){
    
    chr <- combis.df$chr[[i]]
    rt.calc <- combis.df$rt.calc[[i]]
    
    vals <- rt[rt$chroms == chr, rt.calc]
    save(vals, file=paste0(out.dir, "/", chr, "_", rt.type, "_", rt.calc, "_", out.id, ".RData"))
    
    rm(i, chr, rt.calc, vals)
    
  }
  
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()