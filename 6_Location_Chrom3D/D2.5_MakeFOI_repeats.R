################################################################################
# Make foi list of repeats compatible for calculating densities at radial windows
# Mac, R/3.5.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D/out_Chrom3DVsRepeats"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = wk.dir
### OTHER SETTINGS #############################################################
featurefile = paste0(wk.dir, "/Repeat_rankingbyAge/plot_GiorPubl372rankrepFamilies.csv")
out.name = "GiorPubl"
maxgroup = 20L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
feature <- fread(file=featurefile, header=TRUE, data.table=FALSE, sep=",",
                 stringsAsFactors=FALSE)[["repName"]]

counter <- seq(from=1, to=maxgroup, by=1)
counter.len <- length(counter)

# Young repeats
foi.list <- sapply(X=counter, simplify=FALSE, FUN=function(count){
  tail(x=feature, n = count)
})
names(foi.list) <- as.character(counter)
save(foi.list, file=paste0(out.dir, "/", out.name, "_youngRep", 
                           counter.len, ".RData"))

# Old repeats
foi.list <- sapply(X=counter, simplify=FALSE, FUN=function(count){
  head(x=feature, n = count)
})
names(foi.list) <- as.character(counter)

save(foi.list, file=paste0(out.dir, "/", out.name, "_oldRep", 
                           counter.len, ".RData"))

# rm(list=ls())

