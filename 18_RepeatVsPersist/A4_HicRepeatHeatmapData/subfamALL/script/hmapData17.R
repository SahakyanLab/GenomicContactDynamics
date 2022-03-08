################################################################################
# Generate PREELMTISSDYN.MX per chr containing the number of contacts with non-0
# shared number of repeats (or minimum repeat count of the pair of contacting
# regions, Fr^rij) 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/19_RepeatVsPersist")
  } else if (whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")

rep.group = "subfamALL" # "fam" | "subfam" | subfam6
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
minelm.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group)
out.dir = paste0(wk.dir, "/out_HicRepeatHeatmapData/", rep.group)
### OTHER SETTINGS #############################################################
# Age rank identifier
out.name = "subfamALL"
gcb = "min2Mb"
chr = "chr17" #paste("chr", c(1:22, "X"), sep="")
# Max number is number of repeats in age rank
nCPU = 1L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
library(compiler)
source(paste0(lib, "/UTL_doPar.R"))
## FUNCTION ####################################################################
HicRepeatHeatmapData <- function(
  # Vector of elements in agerank
  agerank = agerank,
  minelm.dir = paste0(wk.dir, "/out_HicRepeatExploration"),
  out.dir = paste0(wk.dir, "/out_HicRepeatHeatmapData"),
  gcb = "min2Mb",
  chr = "chr1", #paste("chr", c(1:22, "X"), sep="")
  # Max = number of repeat elements
  nCPU = 30L,
  # Age rank identifier
  out.name = "GiorPubl"
  #combineOut = FALSE
){
  
  agerank <-  as.character(agerank)
  agerank.len <- length(agerank)
  #if(combineOut==TRUE){
  # Initialize final matrix following ELMTISSDYN.MX format
  #  num.non0cont.MX <- matrix(data=0L, nrow=agerank.len + 1L, ncol=21)
  #  dimnames(num.non0cont.MX) <- list(c("num.contact", agerank), 1:21) 
  #}
  load(paste0(minelm.dir,"/", chr, "_MinElm_", gcb, ".RData"))
  num.contacts <- sapply(X=1:21, FUN=function(ntis){
    test2 <- (MINELM.MX[,"ntis"]==ntis)
    sum(test2)
  })
  # Agerank repeat missing in MINELM.MX of this chr (this varies among chr)
  el.missing.minelm <- agerank[!agerank%in%colnames(MINELM.MX[,-1])]
  
  toExport <- c("agerank", "gcb", "MINELM.MX", "el.missing.minelm")
  #### PARALLEL EXECUTION #########
  PREELMTISSDYN.MX <- foreach(itr=isplitVector(1:agerank.len, chunks=nCPU), 
                              .combine="rbind", .inorder=TRUE, 
                              .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      element <- agerank[i]
      print(element, quote=FALSE)
      
      if(element%in%el.missing.minelm){
        non0cont.each.ntis <- rep(0, times=21)
        print("All 0s.")
        return(non0cont.each.ntis)
      } else {
        test1 <- (MINELM.MX[,element]>=1)
        non0cont.each.ntis <- rep(NA,times=21)
        for(ntis in 1:21){
          test2 <- (MINELM.MX[,"ntis"]==ntis)
          # Number of contacts with non0 mininum repeat count
          non0cont.each.ntis[ntis] <- sum(test2 & test1)
          # Sum of minimum repeat counts 
          #non0cont.each.ntis[ntis] <- sum(MINELM.MX[test2, element],)
          rm(test2)
        }#; rm(test1)
        return(non0cont.each.ntis)
      }
    })
    return(do.call("rbind", chunk))
  }
  ### END OF PARALLEL EXECUTION ###
  
  rm("MINELM.MX"); gc()
  
  PREELMTISSDYN.MX <- rbind(num.contacts, PREELMTISSDYN.MX)
  rm(num.contacts)
  if(nrow(PREELMTISSDYN.MX)==(agerank.len+1L) & ncol(PREELMTISSDYN.MX)==21L){
    dimnames(PREELMTISSDYN.MX) <- list(c("num.contacts", agerank), 1:21)
  } else {
    stop("Checkpoint 1.")
  }
  save(x=PREELMTISSDYN.MX, file=paste0(out.dir, "/", chr, "_", gcb, "_",
                                       out.name, "_PreElmTissDyn.RData"))
  rm("PREELMTISSDYN.MX"); gc()
  print(paste0(chr, ":Heatmap preliminary data generated!"))
}
################################################################################
HicRepeatHeatmapData <- cmpfun(HicRepeatHeatmapData, options=list(suppressUndefined=TRUE))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.dir <- paste0(out.dir, "/", rep.group)
if( !dir.exists(out.dir) ){
  dir.create(out.dir)
}

col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                    header=TRUE, stringsAsFactors=FALSE)[,col.nme]
agerank <- as.character(agerank)

#for(chr in chr.v){
  
  print(paste0(gcb, "_", rep.group, "_", chr, "..."), quote=FALSE)

  HicRepeatHeatmapData(
    agerank=agerank,
    minelm.dir=minelm.dir,
    out.dir=out.dir,
    gcb=gcb,
    chr=chr, 
    nCPU=nCPU,
    out.name=out.name
  )
  
#}

# rm(list=ls()); gc()
