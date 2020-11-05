################################################################################
# HicRepeatExploration 
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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "fam" # "fam" | "subfam"
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
binRep.dir = paste0(wk.dir, "/out_RepeatOverlapPerBin/", rep.group)
out.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group)
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = "chr7" 
nCPU = 4 # Number of contacts
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(itertools)
library(doParallel)
library(compiler)
source(paste0(lib, "/UTL_doPar.R"))
### FUNCTION ###################################################################
HicRepeatExploration <- function(
  # Number of CPU cores to be used
  nCPU = 4, # optimal for the task on bucephalus (takes ~5GB per core for chr1)
  # Chromosome for which the data are to be pulled
  chr = 21,
  # PERSIST.MX directory 
  PersistPath = persist.dir,
  # BINREP.MX directory
  BinRepPath = binRep.dir, 
  # directory for intermediate and final files
  out.dir = out.dir,
        #"/Volumes/Data/Database/HiC_features_GSE87112_RAWpc",
  # Feature database filename suffix and will be used for saved outputs
  suffix = "min2Mb"
){
  
id <- paste0(chr,"_",suffix)

load(paste0(binRep.dir, "/", suffix,"_", chr, "_BinRep.RData"))
load(paste0(persist.dir, "/", chr, "_Persist_", suffix,".RData"))

############################
element.names <- names(BINREP.MX[1,-(1:3)])
hits.len <- length(PERSIST.MX$hits[,1])
ij.mx <- cbind(i=PERSIST.MX$hits[,"i"], j=PERSIST.MX$hits[,"j"])
ntis.vec <- PERSIST.MX$ntis

rm(PERSIST.MX)
gc()

elem.names.counter <- 0L

###---------------------------
for(element in element.names){
  
  print(element, quote=FALSE)
  
  #### FOREACH EXECUTION #########
  element.count <- foreach(itr=isplitVector(1:hits.len, chunks=nCPU),
                           .combine="c", .inorder=TRUE,
                           .export=c("BINREP.MX","ij.mx","element"),
                           .noexport=ls()[!ls()%in%c("BINREP.MX","ij.mx","element")]
  ) %op% {
    
    element.count.chunk <- sapply(itr,
                                  FUN=function(itr){
                                    min(
                                      BINREP.MX[which(BINREP.MX[,"bins"] %in%
                                                        ij.mx[itr,]),element]
                                    )
                                  },
                                  simplify=TRUE, USE.NAMES=FALSE)
    return(element.count.chunk)
  }
  ### END OF FOREACH EXECUTION ###
  
  # Dumped per element, then collected later, for memory efficiency.
  
  # used a counter for saving because some repeat subfamilies have characters 
  # inappropriate for saving
  elem.names.counter <- elem.names.counter + 1
  save(element.count, file=paste0(out.dir,"/",chr, "_MinElmSOURCE_",
                                  suffix,"_",as.character(elem.names.counter),".RData"))
  
  rm(element.count); gc()
  
}
###---------------------------

rm(ij.mx); gc()

# Saving all data for min(repeat counts) in each interacting pairs of chr loci
MINELM.MX <- matrix( NA, nrow=hits.len, ncol=length(element.names)+1 )
dimnames(MINELM.MX)[[2]] <- c("ntis", element.names)
MINELM.MX[,"ntis"] <- ntis.vec; rm(ntis.vec)
elem.names.len <- length(element.names)
for(elem.names.counter in 1:elem.names.len){
  SRC<-paste0(out.dir,"/",chr,"_MinElmSOURCE_",suffix,
              "_",as.character(elem.names.counter),".RData")
  load(SRC); file.remove(SRC)
  MINELM.MX[,element.names[elem.names.counter]] <- element.count
  rm(element.count)
  gc()
}
save(MINELM.MX, file=paste0(out.dir,"/",chr,"_MinElm_",suffix,".RData"))

############################

print("HicRepeatExploration is DONE!", quote=FALSE)

}
################################################################################
HicRepeatExploration <- cmpfun(HicRepeatExploration, options=list(suppressUndefined=TRUE))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, "..."), quote=FALSE)

for(chr in chr.v){
  
  HicRepeatExploration(
    # Number of CPU cores to be used
    nCPU = nCPU,
    # Chromosome for which the data are to be pulled
    chr = chr,
    # PERSIST.MX directory 
    PersistPath = persist.dir,
    # BINREP.MX directory
    BinRepPath = binRep.dir, 
    # directory for intermediate and final files of HiCRepeatExploration
    out.dir = out.dir,
    # Feature database filename suffix
    suffix = gcb
  )
  
}
################################################################################
#changes made with Alex's original script
##deleted initialExplPlots section
##removed initialExplPlots = TRUE and mobDNAExplPlots = TRUE arguments
##element.names <- names(BINREP.MX[1,-(1:3)])
##load(paste0(out.dir,"/",chr,"_BinRep_Subfam",suffix,".RData"))
##skipped generation of plots, too much for repNames

