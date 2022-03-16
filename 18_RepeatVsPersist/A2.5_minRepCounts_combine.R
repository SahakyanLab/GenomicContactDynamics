################################################################################
# Boxplot of minimum repeat count per Cp combining all chromosome data
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")

# Metric
metric = "_skewrep" # _skewrep | "" for minrep

rep.group = "subfamALL" #"subfamALL" # "fam" | "subfam"
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
out.dir = paste0(wk.dir, "/out_minRepCounts/", rep.group, metric)
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
ntis.v = 1:21
# Number of repeat elements (372/56/62)
nCPU = 2L # 15G each --> set to 50
# Age rank identifier
out.name = paste0(rep.group, metric) #"subfamALL_minRepCounts" #"GiorPubl_minRepCounts"
plotOnly = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/compareTwoDist.R"))
source(paste0(wk.dir, "/lib/makeMinRepPlot.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, " ", rep.group, "..."), quote=FALSE)

if(plotOnly==FALSE){
  
  col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
  agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                      header=TRUE, stringsAsFactors=FALSE)[,col.nme]; rm(col.nme)
  agerank <- as.character(agerank)
 
  # Initialize chrALL list
  lst <- vector("list", length(ntis.v))
  names(lst) <- as.character(ntis.v)
  MINREPCOUNTS.ALL <- rep(list(lst), length(agerank))
  names(MINREPCOUNTS.ALL) <- agerank
  
  for(chr in chr.v){
    
    load(paste0(out.dir, "/", chr, "_", gcb, "_", out.name, ".RData"))
    
    if( !identical(agerank, names(MINREPCOUNTS)) ){
      
      rm(MINREPCOUNTS)
      stop(paste0(chr, ": Faulty MINREPCOUNTS."))
      
    } else {
      
      elements <- agerank
      elements.len <- length(elements)
      
    }
    
    toExport <- c("elements", "MINREPCOUNTS", "MINREPCOUNTS.ALL")
    
    # REMOVE
    print("#### PARALLEL EXECUTION #########", quote=FALSE)
    
    #### PARALLEL EXECUTION #########
    
    # List of lists
    lst.temp <- foreach( itr=isplitVector(1:elements.len, chunks=nCPU), 
                         .inorder=TRUE, .combine="c", .export=toExport, 
                         .noexport=ls()[!ls()%in%toExport]
                                       
    ) %op% {
      
      chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
        
        elm <- elements[i]
        #print(elm)
        
        lst.ntis <- sapply(X=ntis.v, simplify=FALSE, FUN=function(ntis){
          
          v <- c(MINREPCOUNTS.ALL[[elm]][[ntis]], MINREPCOUNTS[[elm]][[ntis]])
          
          if( is.null(v) | length(v)==0 ){
            
            return(NULL)
           
          } else {
            
            agg <- aggregate(x=v, by=list(names(v)), FUN=sum)
            return(setNames(object=agg$x, nm=agg$Group.1))
            
          }
          
        })
        
        names(lst.ntis) <- ntis.v
        
        return(lst.ntis)
        
      }) # itr sapply end
      return(chunk)
      
    }
    
    ### END OF PARALLEL EXECUTION ###
    
    # elements == agerank
    if( length(lst.temp) != length(elements) ){
      rm(MINREPCOUNTS.ALL)
      stop("Lengths not equal.")
    } else {
      MINREPCOUNTS.ALL <- setNames(object=lst.temp, nm=elements)
    }
    
    rm(MINREPCOUNTS, elements, lst.temp); gc()
    print(paste0(chr, " done!"), quote=FALSE)
    
  } # chr.v for loop end
  
  MINREPCOUNTS <- MINREPCOUNTS.ALL
  rm(MINREPCOUNTS.ALL)
  
  save(MINREPCOUNTS, 
       file=paste0(out.dir, "/chrALL_", gcb, "_", out.name, ".RData"))
  
} else {
  
  # Load MINREPCOUNTS
  load(file=paste0(out.dir, "/chrALL_", gcb, "_", out.name, ".RData"))

}

makeMinRepPlot(MINREPCOUNTS, ntis.v, affix=paste0("chrALL_", gcb, "_", out.name), out.dir, nCPU=nCPU)

# rm(list=ls()); gc()