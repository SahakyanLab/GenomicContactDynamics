################################################################################
# 
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

metric = "minrep"
src.dir = paste0(wk.dir, "/z_ignore_git/out_minRepCounts/subfamALL_", metric, "_atleast2sumrep")
out.dir = paste0(wk.dir, "/out_combine") #paste0(wk.dir, "/out_metric_estimate")
### OTHER SETTINGS #############################################################
src.nme = paste0("chrALL_min2Mb_subfamALL_", metric, "Counts")
Cps.forCalc = 19:21 # Cp values to consider when calculating median non-zero count fraction
out.id = paste0("estimate_Cp", Cps.forCalc[[1]], "To", tail(Cps.forCalc, n=1))
# GROUP.CLASS <- list(
#   not.transposon = c("Low_complexity", "RNA", "rRNA", "Satellite", "scRNA", 
#                      "Simple_repeat", "snRNA", "srpRNA", "tRNA"),
#   DNA.transposon = c("DNA", "DNA?"), 
#   retro.transposon = c("LINE", "SINE", "LTR", "RC", "LINE?", "SINE?", "LTR?"),
#   not.classified = c("Other", "Unknown", "Unknown?")
# )
# group.cols = setNames(object=c("#7E6148FF", "#3C5488FF", "#00A087FF", "black"), 
#                       nm=names(GROUP.CLASS))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(src.dir, "/", src.nme, ".RData"))
elements <- names(MINREPCOUNTS)

# Deal with NULLS, these are for Cp values without contacts with at least 2 sites
# of the element; convert nulls to c(`0`=NA)
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

# Metric estimate = median of values from Cps.forCalc

CP.ESTIMATE <- MINREPCOUNTS
for(elm in elements){
  
  cpvals <- MINREPCOUNTS[[elm]][ c(as.character(Cps.forCalc)) ]
  cpvals <- lapply(cpvals, FUN=stack)
  cpvals <- do.call("rbind", cpvals)
  setnames(cpvals, old=c("values", "ind"), new=c("count", "metric.val"))
  cpvals <- cpvals[!is.na(cpvals$count),]
  
  if( nrow(cpvals) > 0 ){
    cpvals <- Map(f=rep, as.numeric(as.character(cpvals$metric.val)), times=cpvals$count)
    cpvals <- median(unlist(cpvals))
  } else {
    cpvals <- NA
  }
  
  CP.ESTIMATE[[elm]] <- cpvals
    
}

## PLOT data

df <- stack(unlist(CP.ESTIMATE))
df$ind <- as.character(df$ind)
setnames(df, old="ind", new="repName")

save(df, file=paste0(out.dir, "/", metric, "_", out.id, ".RData"))

# rm(list=ls()); gc()