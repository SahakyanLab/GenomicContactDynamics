################################################################################
# Per subfamily, determine number of contacts of given Cp or Cp range with at 
# least 1 shared number of site. Cps.forCalc is a list of Cp or Cp range to use.
# If Cps.forCalc pertains to more than 1 Cp value, the number would be determined
# using all contacts within Cp range.
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

src.dir = paste0(wk.dir, "/z_ignore_git/out_minRepCounts/subfamALL_minrep") #, "_atleast2sumrep")
out.dir = paste0(wk.dir, "/z_ignore_git/out_combine") 
### OTHER SETTINGS #############################################################
src.nme = "chrALL_min2Mb_subfamALL_minrepCounts"
Cps.forCalc = list(dyn1T3=1:3, per19To21=19:21) # Cp values to consider when calculating median non-zero count fraction
out.id = paste0("metric_nonzero_Cp", names(Cps.forCalc)[[1]], "_", names(Cps.forCalc)[[2]])
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(src.dir, "/", src.nme, ".RData"))
elements <- sort(names(MINREPCOUNTS), decreasing=F)

# # Deal with NULLS, these are for Cp values without contacts with at least 2 sites
# # of the element; convert nulls to c(`0`=NA)
# for(elm in elements){
#   
#   MINREPCOUNTS[[elm]] <- sapply(names(MINREPCOUNTS[[elm]]), simplify=F, FUN=function(cp.nme){
#     counts <- MINREPCOUNTS[[elm]][[cp.nme]]
#     if( is.null(counts) ){
#       message(elm)
#       return(c(`0`=NA)) 
#     } else {
#       return(counts)
#     }
#   })
#   
# }

if( unique(lengths(MINREPCOUNTS)) != 21 ){
  stop("Missing Cps.")
}

#

Cps.forCalc.len <- length(Cps.forCalc)
NONZERO.COUNTS <- list()

for( Cps.i in 1:Cps.forCalc.len){

  # # Sum of contact counts combining Cps.forCalc
  # 
  # SUM.COUNTS <- MINREPCOUNTS
  # for(elm in elements){
  #   SUM.COUNTS[[elm]] <- sum(unlist( MINREPCOUNTS[[elm]][ Cps.forCalc[[Cps.i]] ] ), na.rm=T)
  #   #SUM.COUNTS[[elm]] <- lapply(MINREPCOUNTS[[elm]], FUN=sum, na.rm=T)
  # }
  # 
  # if( any( !is.finite(unlist(SUM.COUNTS)) ) ){
  #   stop("Non-finite sum.")
  # }
  
  # Contact count with >=1 shared site combining Cps.forCalc
  
  NONZERO.COUNTS.Cps.i <- MINREPCOUNTS
  for(elm in elements){
    
    nonzero.counts <- lapply(MINREPCOUNTS[[elm]][ Cps.forCalc[[Cps.i]] ], FUN=function(cp.counts){
      cp.counts <- cp.counts[names(cp.counts) != "0"]
      return(cp.counts)
    })
    NONZERO.COUNTS.Cps.i[[elm]] <- sum(unlist( nonzero.counts ), na.rm=T)
    
  }
  
  NONZERO.COUNTS.Cps.i <- stack(NONZERO.COUNTS.Cps.i)
  NONZERO.COUNTS[[Cps.i]] <- NONZERO.COUNTS.Cps.i$values

} # Cps.forCalc.len for loop end

if( !is.null(names(Cps.forCalc)) ){
  names(NONZERO.COUNTS) <- names(Cps.forCalc)
}

NONZERO.COUNTS <- do.call("cbind", NONZERO.COUNTS)
dimnames(NONZERO.COUNTS)[[1]] <- elements

mx <- NONZERO.COUNTS

save(mx, file=paste0(out.dir, "/", out.id, ".RData"))

# rm(list=ls()); gc()