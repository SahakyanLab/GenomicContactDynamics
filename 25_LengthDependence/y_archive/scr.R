################################################################################
# Check influence of gap between contacting bins to the hybridisation trend 
# deva, R/3.6.0-newgcc, gcc/4.9.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/16_LengthDependence"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/16_LengthDependence"
    data.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = hyb.dir = data.dir 
out.dir = paste0(wk.dir, "/out")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1"
kmer.len = 7
# In terms of bin
splitLength = 550
nCPU = 1L
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
library(RColorBrewer)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/HiCHybridPlot.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# PERSIST.MX
load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))

# HYB.MX
load(file=paste0(data.dir, "/", chr, "_Hyb", kmer.len, "_", gcb, ".RData"))
gap.bin <- PERSIST.MX$hits$j-PERSIST.MX$hits$i

# Get Cp data
cp.v <- PERSIST.MX$ntis
rm(PERSIST.MX); gc()

# Numbers
tot.ij <- length(HYB.MX[1,])
nonNA.ij.TF <- !is.na(HYB.MX[1,])

out.name <- paste0(chr, "_" , gcb, "_kmer", kmer.len, "_splitLength", splitLength)

toExport <- c("tot.ij", "cp.v", "out.dir", "out.name", "HYB.MX", "nonNA.ij.TF",
              "gap.bin", "splitLength")

foreach(grp=c("short", "long"), .inorder=FALSE, 
        .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  if(grp=="short"){
    TF <- gap.bin <= splitLength
  } else if(grp=="long"){
    TF <- gap.bin > splitLength
  } 
  
  incl.ij.TF <- nonNA.ij.TF & TF
  incl.ij.Perc <- round(x=sum(incl.ij.TF)/tot.ij*100, digits=2)
  percPerCp <- table(cp.v[incl.ij.TF])/tot.ij*100
  percPerCp <- round(x=percPerCp[as.character(1:21)], digits=4)
  percPerCp <- paste(percPerCp, collapse=" ")
  
  HiCHybridPlot(
    out.dir=out.dir,
    out.name=paste0(out.name, "_", incl.ij.Perc, "Perc", grp),
    plottitle=paste0(out.name, "_", tot.ij, "totij_", 
                     round(100*sum(nonNA.ij.TF)/tot.ij, digits=2), 
                     "%nonNA_", incl.ij.Perc, "%", grp, "\n", percPerCp),
    HYB.MX=HYB.MX[,incl.ij.TF],
    # Corresponding order to HYB.MX (column)
    cp.v=cp.v[incl.ij.TF],
    # Corresponding order to HYB.MX (row)
    label.v=c("T.H.E. Gfree, kcal/mol",
              "T.H.E. sd(discordance)",
              "T.H.E. sum(abs(discordance))")
  )
  
  print(paste0(grp, " done!"), quote=FALSE)
  
}
# rm(list=ls())
  