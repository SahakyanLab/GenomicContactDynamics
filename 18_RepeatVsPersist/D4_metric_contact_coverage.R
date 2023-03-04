################################################################################
# Calculate per subfamily the percentage of contacts with measured metrics so
# contacts with at least 2 sites relative to number of long-range contacts per
# Cp. 
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
out.dir = paste0(wk.dir, "/out_combine") #paste0(wk.dir, "/out_metric_contact_coverage")
ijcount.file = paste0(wk.dir, "/min2Mb_ij_2Mb_Cp1To21_ijcount.txt")
### OTHER SETTINGS #############################################################
Cps = 1:21 # Same order as ijcount.file (top to bottom)
out.id = "contact_coverage"
src.nme = paste0("chrALL_min2Mb_subfamALL_", metric, "Counts")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
ij.count <- setNames(as.numeric(readLines(ijcount.file)), nm=Cps)

Cps.len <- length(Cps)

load(paste0(src.dir, "/", src.nme, ".RData"))
elements <- sort(names(MINREPCOUNTS), decreasing=F)
elements.len <- length(elements)

# Deal with NULLS, these are for Cp values without contacts with at least 2 sites
# of the element; convert nulls to c(`0`=NA)

for(elm in elements){
  
  MINREPCOUNTS[[elm]] <- lapply(MINREPCOUNTS[[elm]], FUN=function(counts){
    if(is.null(counts)){
      return(c(`0`=NA)) 
    } else {
      return(counts)
    }
  })
  
}

# Sum of contact counts across values per Cp per element

SUM.COUNTS <- MINREPCOUNTS
for(elm in elements){
  SUM.COUNTS[[elm]] <- lapply(MINREPCOUNTS[[elm]], FUN=sum)
  SUM.COUNTS[[elm]] <- do.call("cbind", SUM.COUNTS[[elm]])
}

SUM.COUNTS <- do.call("rbind", SUM.COUNTS)
# if( identical(dimnames(SUM.COUNTS)[[2]], as.character(Cps)) ){
#   # Need to transpose cause R stores its matrix information in column major order
#   # https://stackoverflow.com/questions/20596433/how-to-divide-each-row-of-a-matrix-by-elements-of-a-vector-in-r
#   SUM.COUNTS <- t( t(SUM.COUNTS) / as.numeric(ij.count) * 100 )
# }

dimnames(SUM.COUNTS)[[1]] <- elements

mx <- SUM.COUNTS
save(mx, file=paste0(out.dir, "/", metric, "_", out.id, ".RData"))


# rm(list=ls()); gc()