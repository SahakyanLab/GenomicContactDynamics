################################################################################
# Using BINREP.MX, get shared number of each repeat per contact, per chr.
# The shared number is the minimum/smaller repeat count between the two
# contacting regions.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/18_RepeatVsPersist")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" 
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")

rep.group = "subfam" # "fam" | "subfam"
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
binRep.dir = paste0(wk.dir, "/out_RepeatOverlapPerBin/", rep.group)
out.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group, "ALL")
source.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group, "ALL_SOURCE")
elementlistPath = paste0(wk.dir, "/out_makeElementsList/", rep.group)
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr7" 
nCPU = 1L # Number of contacts
# Index of element in elementlistPath; has to be an integer so add 'L', NA if
# not element-wise running
element.ind = NA #elementREPLACE
makeMinElmSOURCE = FALSE
makeMINELMMX = TRUE
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

  PersistPath = 'PERSIST.MX directory',
  BinRepPath = 'BINREP.MX directory', 
  out.dir = 'Directory for intermediate and final files',
  elementlistPath = 'Path to txt file of unique element names in BINREP.MX',
  element.ind = 'If not NA, run code only for this element. If NA, run code
                 for all element names in BINREP.MX',
  nCPU = 'Number of cores for parallel execution; depends on number of contacts', 
  chr = 'chromosome',
  suffix = 'contact gap id i.e. min2Mb or min05Mb',
  makeMinElmSOURCE = TRUE,
  makeMINELMMX = FALSE
  
){
  
id <- paste0(chr,"_",suffix)

load(paste0(binRep.dir, "/", chr, "_BinRep_", suffix, ".RData"))
#load(paste0(binRep.dir, "/", chr, "_BinRep.RData"))
load(paste0(persist.dir, "/", chr, "_Persist_", suffix,".RData"))

############################
if( is.na(element.ind) ){
  
  element.names <- names(BINREP.MX[1,-(1:3)])
  elem.names.counter <- 0L
  
} else if( is.integer(element.ind) & length(element.ind)==1 ){
  
  elementlistPath = paste0(elementlistPath, "/", suffix, "_", chr, "_elements.txt")
  tmp <- readLines(con=elementlistPath) 
  
  if( element.ind>length(tmp) ){
    stop("HicRepeatExploration(): element.ind out of bounds.")
  } else {
    element.names <- tmp[[element.ind]]
  }
  rm(tmp)
  
  elem.names.counter <- element.ind-1L
  
} else {
  stop("HicRepeatExploration(): Invalid argument.")
}

hits.len <- length(PERSIST.MX$hits[,1])
ij.mx <- cbind(i=PERSIST.MX$hits[,"i"], j=PERSIST.MX$hits[,"j"])
ntis.vec <- PERSIST.MX$ntis

rm(PERSIST.MX)
gc()

###---------------------------
if(makeMinElmSOURCE){
  
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
    save(element.count, file=paste0(source.dir, "/", chr, "_MinElmSOURCE_", suffix,
                                    "_", as.character(elem.names.counter), ".RData"))
    
    rm(element.count); gc()
    
  }

}
###---------------------------
rm(ij.mx)
gc()

if(makeMINELMMX){
  
  # Saving all data for min(repeat counts) in each interacting pairs of chr loci
  MINELM.MX <- matrix( NA, nrow=hits.len, ncol=length(element.names)+1 )
  dimnames(MINELM.MX)[[2]] <- c("ntis", element.names)
  MINELM.MX[,"ntis"] <- ntis.vec
  rm(ntis.vec)
  elem.names.len <- length(element.names)
  
  for(elem.names.counter in 1:elem.names.len){
    
    print(elem.names.counter, quote=FALSE)
    SRC <- paste0(source.dir, "/", chr, "_MinElmSOURCE_", suffix, "_", 
                  as.character(elem.names.counter),".RData")
    load(SRC)
    #file.remove(SRC)
    MINELM.MX[,element.names[elem.names.counter]] <- element.count
    rm(element.count)
    gc()
    
  }
  save(MINELM.MX, file=paste0(out.dir, "/", chr, "_MinElm_", suffix, ".RData"))
  
}

############################

print("HicRepeatExploration is DONE!", quote=FALSE)

}
################################################################################
HicRepeatExploration <- cmpfun(HicRepeatExploration, options=list(suppressUndefined=TRUE))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, "..."), quote=FALSE)

HicRepeatExploration(
  
  PersistPath=persist.dir,
  BinRepPath=binRep.dir, 
  out.dir=out.dir,
  elementlistPath=elementlistPath,
  element.ind=element.ind,
  nCPU=nCPU,
  chr=chr,
  suffix=gcb,
  makeMinElmSOURCE=makeMinElmSOURCE,
  makeMINELMMX=makeMINELMMX
  
)
################################################################################
#changes made with Alex's original script
##deleted initialExplPlots section
##removed initialExplPlots = TRUE and mobDNAExplPlots = TRUE arguments
##element.names <- names(BINREP.MX[1,-(1:3)])
##load(paste0(out.dir,"/",chr,"_BinRep_Subfam",suffix,".RData"))
##skipped generation of plots, too much for repNames

# rm(list=ls()); gc()