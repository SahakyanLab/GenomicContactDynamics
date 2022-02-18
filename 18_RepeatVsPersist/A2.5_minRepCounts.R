################################################################################
# Boxplot of minimum repeat count per Cp 
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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/18_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")

rep.group = "subfamALL" #"subfamALL" # "fam" | "subfam"
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
minelm.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group)
out.dir = paste0(wk.dir, "/out_minRepCounts/", rep.group)
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "CHRREPLACE"
ntis.v = 1:21
# Number of repeat elements (372/56/62)
nCPU = 1L # 15G each --> set to 50
# Age rank identifier
out.name = "subfamALL_minRepCounts" #"subfamALL_minRepCounts" #"GiorPubl_minRepCounts"
#plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(wk.dir, "/lib/makeMinRepPlot.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( !dir.exists(out.dir) ){
  dir.create(out.dir)
}

print(paste0(gcb, " ", rep.group, "..."), quote=FALSE)

#if(plotOnly==FALSE){
  
  col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
  agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                      header=TRUE, stringsAsFactors=FALSE)[,col.nme]; rm(col.nme)
  agerank.len <- length(agerank)

  #chr.v.len <- length(chr.v)
  #for(ch in 1:chr.v.len){
    
    #chr <- chr.v[ch]
    load(paste0(minelm.dir,"/", chr, "_MinElm_", gcb, ".RData"))
    print(paste0(chr, " MINELM.MX loaded."), quote=FALSE)
    
    elements <- intersect(agerank, dimnames(MINELM.MX)[[2]][-1])
    elements.len <- length(elements)
    
    MINELM.MX <- MINELM.MX[,c("ntis", elements)]
    
    #if(ch==1){
      
      # Initialize list
      #lst <- vector("list", length(ntis.v))
      #names(lst) <- as.character(ntis.v)
      #MINREPCOUNTS <- rep(list(lst), agerank.len)
      MINREPCOUNTS <- vector("list", agerank.len)
      names(MINREPCOUNTS) <- agerank
      
    #}
    
    toExport <- c("MINELM.MX", "elements", "MINREPCOUNTS")
    
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
          
          test <- MINELM.MX[,"ntis"]==ntis
          
          if(sum(test)>0){
            
            #v <- c(MINREPCOUNTS[[elm]][[ntis]], table(MINELM.MX[test,elm]))
            #agg <- aggregate(x=v, by=list(names(v)), FUN=sum)
            #return(setNames(object=agg$x, nm=agg$Group.1))
            return( table(MINELM.MX[test,elm]) )
            
          } else {
            
            warning(paste0(elm, ": No Cp=", ntis, " in ", chr, "."))
            #return( MINREPCOUNTS[[elm]][[ntis]] )
            return(NULL)
     
          }
          
        })
        
        names(lst.ntis) <- ntis.v
        
        return(lst.ntis)
        
      }) # itr sapply end
      return(chunk)
      
    }
    
    ### END OF PARALLEL EXECUTION ###
    
    MINREPCOUNTS[elements] <- lst.temp
    rm(lst.temp); gc()
    
    rm(MINELM.MX, elements); gc()
      
  #} # chr.v for loop end
  
  save(MINREPCOUNTS, 
       file=paste0(out.dir, "/", chr, "_", gcb, "_", out.name, ".RData"))
  
#} else {
  
#  # Load MINREPCOUNTS
#  load(file=paste0(out.dir, "/", chr, "_", gcb, "_", out.name, ".RData"))

#}

#makeMinRepPlot(MINREPCOUNTS, ntis.v, affix=paste0(chr, "_", gcb, "_", out.name), out.dir)

# rm(list=ls()); gc()