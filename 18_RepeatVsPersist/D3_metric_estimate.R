################################################################################
# Per subfamily, determine median or mean (estimate.fnx argument) value of 
# metrics of contacts of given Cp or Cp range with at least 1 shared number of 
# site. Cps.forCalc is a list of Cp or Cp range to use. If Cps.forCalc pertains 
# to more than 1 Cp value, the number would be determined using all contacts 
# within Cp range. Compare distributions of each Cp or Cp range in Cps.forCalc. 
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")

metric = "skewrep"
src.dir = paste0(wk.dir, "/z_ignore_git/out_minRepCounts/subfamALL_", metric) #, "_atleast2sumrep")
out.dir = paste0(wk.dir, "/z_ignore_git/out_combine") #paste0(wk.dir, "/out_metric_estimate")
### OTHER SETTINGS #############################################################
src.nme = paste0("chrALL_min2Mb_subfamALL_", metric, "Counts")
Cps.forCalc = list(dyn1T3=1:3, per19To21=19:21) # Cp values to consider when calculating median non-zero count fraction
estimate.fnx = "median" #"mean" | "median" 
out.id = paste0("metric_estimate", estimate.fnx, "_Cp", names(Cps.forCalc)[[1]], 
                "_", names(Cps.forCalc)[[2]])
out.dir = paste0(wk.dir, "/z_ignore_git/out_combine") #/estimates_", estimate.fnx)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(foreach)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(src.dir, "/", src.nme, ".RData"))
#MINREPCOUNTS <- MINREPCOUNTS[1:5] # REMOVE
elements <- names(MINREPCOUNTS)

# Confirm that only metric values that exist in at least one contact are recorded
# in MINREPCOUNTS

counts <- unlist(MINREPCOUNTS)
if( any(counts <= 0) ){
  rm(MINREPCOUNTS)
  stop("<= 0 counts")
}
rm(counts)
  
# Deal with NULLS, these are for Cp values without contacts with at least 2 sites
# of the element; convert nulls to c(`0`=NA)
# c(`0`=NA) is converted to 0 metric val with 0 counts and are removed
for(elm in elements){

  MINREPCOUNTS[[elm]] <- sapply(names(MINREPCOUNTS[[elm]]), simplify=F, FUN=function(cp.nme){
    counts <- MINREPCOUNTS[[elm]][[cp.nme]]
    if(is.null(counts)){
      return(c(`0`=NA))
    } else {
      return(counts)
    }
  })

}

if( unique(lengths(MINREPCOUNTS)) != 21 ){
  stop("Missing Cps.")
}

# Metric estimate = median/mean of values from Cps.forCalc

#mx <- matrix(data=NA, nrow=length(elements), ncol=length(Cps.forCalc) + 1,
#             dimnames=list(elements, c(names(Cps.forCalc), "pval")))

Cps.forCalc.len <- length(Cps.forCalc)

ESTIMATE <- list()
for(elm in elements){
  
  # VALCOUNTS is a list containing dataframe of metric values and counts per Cps.forCalc, 
  VALCOUNTS <- foreach(Cps.i=1:Cps.forCalc.len, .combine="list", .inorder=T
  
  ) %do% {
    
    cpvals <- MINREPCOUNTS[[elm]][ c(as.character(Cps.forCalc[[Cps.i]])) ]
    # df of values and counts
    cpvals <- do.call("rbind", lapply(cpvals, FUN=stack))
    setnames(cpvals, old=c("values", "ind"), new=c("count", "metric.val"))
    
    # Aggregate to combine counts of same valuess
    cpvals <- aggregate(x=cpvals$count, by=list(cpvals$metric.val), FUN=sum, na.rm=T)
    setnames(cpvals, old=c("Group.1", "x"), new=c("metric.val", "count"))
    cpvals$metric.val <- as.numeric(as.character(cpvals$metric.val))
    cpvals <- cpvals[cpvals$count > 0,]
    
    return(cpvals)
    
  }
  
  ESTIMATE[[elm]] <- lapply(VALCOUNTS, FUN=function(cpvals){
    
    if( nrow(cpvals) == 0 ){
      return(NA)
    } else {
      
      sum.count <- sum(cpvals$count)
      halfsum.val <- sum.count / 2
      
      max.count <- max(cpvals$count)
      val.max.count <- cpvals$metric.val[ cpvals$count == max.count ]
      
      #
      
      if( (max.count > halfsum.val) & (estimate.fnx == "median") ){
        return(val.max.count)
      } else {
        
        message(paste0(elm, " expand..."))
        
        cpvals <- Map(f=rep, as.numeric(as.character(cpvals$metric.val)), times=cpvals$count)
        
        eval(parse(text=paste0(
          'estimate <- ', estimate.fnx, '(unlist(cpvals))'
        )))
        
        return(estimate)
        
      }
      
    }
      
  })
  
  rm(VALCOUNTS)
  
  message(paste0(elm, " done!"))
  
} # elements for loop end

ESTIMATE <- lapply(ESTIMATE, FUN=function(x) do.call("cbind", x))
ESTIMATE <- do.call("rbind", ESTIMATE)

if( !is.null(names(Cps.forCalc)) ){
  dimnames(ESTIMATE)[[2]] <- names(Cps.forCalc)
}

dimnames(ESTIMATE)[[1]] <- elements

mx <- ESTIMATE
save(mx, file=paste0(out.dir, "/", metric, "_", out.id, ".RData"))

# rm(list=ls()); gc()

# #### PARALLEL EXECUTION #########
# foreach(itr=isplitVector(1:cp.v.len, chunks=nCPU), 
#         .inorder=T, .combine="rbind",
#         .export=toExport, .noexport=ls()[!ls()%in%toExport]
#         
# ) %op% {
#   
# }
# ### END OF PARALLEL EXECUTION ###