################################################################################
# Implementation of plotFeatDensPerRadWindow() to plot density of feature at 
# overlapping radial windows of the model. Densities of features are normalized 
# to the amount of DNA at each window.
# deva, R/3.5.0-newgcc, gcc/4.9.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "IMR90_LMNB1_GSE49341_hg19" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
domXYZRFeat.dir = paste0(wk.dir, "/out_DomainVsFeature")
out.dir = paste0(wk.dir, "/out_FeatDensPerRadWindow")
### OTHER SETTINGS #############################################################
ploidy = "haploid"
# Set radius of bin/window
#Currently, the difference between rVal values ~0.005, (numPoints=1000). I can 
# decrease the radius of window to decrease the degree of overlap between windows.
#0.05?
dr = 0.5
hist.breaks = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.6) 
#seq(from=0, to=6, by=0.5)
regenerateCountData = TRUE
regeneratePlotData = TRUE

periphery.lst <- list(1,0, c(1:0))

# Number of features
nCPU = 5L #~41G
#---------------------------------------
# Contact persistence
#---------------------------------------
out.name = paste0(model.id, "_", ploidy)

cp.list = as.list(as.character(c(1:21)))
names(cp.list) = as.character(cp.list)

master.lst = list(
  
  min2Mb_CpFeat = list(
    DOMXYZRFEATfile=paste0(domXYZRFeat.dir, "/", out.name, "_", 
                           "min2Mb_CpFeat_domXYZRFeat.RData"),
    foi.list=cp.list,
    feat.type=expression(bold( "c"["p"] )),
    foiRData=FALSE,
    multiplier=100L,
    # Cannot generate plot for raw counts because it exceeds the limit
    count.v=c("norm")
  )
  #,
  
  #min05Mb_CpFeat = list(
  #  DOMXYZRFEATfile=paste0(domXYZRFeat.dir, "/", out.name, "_", 
  #                         "min05Mb_CpFeat_domXYZRFeat.RData"),
  #  foi.list=cp.list,
  #  feat.type=expression(bold( "c"["p"] )),
  #  foiRData=FALSE,
  #  multiplier=100L,
  #  # Cannot generate plot for raw counts because it exceeds the limit
  #  count.v = c("norm")
  #)
  #,
  
  #---------------------------------------
  # Repeat elements (hg19)
  #---------------------------------------
  #GiorPubl_youngRep20 = list(
  #  DOMXYZRFEATfile=paste0(domXYZRFeat.dir, "/", out.name, "_", 
  #                         "hg19Repeats_domXYZRFeat.RData"),
  #  foi.list=paste0(repfeat.dir, "/GiorPubl_youngRep20.RData"),
  #  feat.type=expression(bold( "N"["young"] )),
  #  foiRData=TRUE,
  #  multiplier=100L,
  #  count.v = c("raw", "norm")
  #),
  
  #GiorPubl_oldRep20 = list(
  #  DOMXYZRFEATfile=paste0(domXYZRFeat.dir, "/", out.name, "_", 
  #                         "hg19Repeats_domXYZRFeat.RData"),
  #  foi.list=paste0(repfeat.dir, "/GiorPubl_oldRep20.RData"),
  #  feat.type=expression(bold( "N"["old"] )),
  #  foiRData=TRUE,
  #  multiplier=100L,
  #  count.v = c("raw", "norm")
  #)
  
)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
library(foreach)
library(itertools)
library(doParallel)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/BINorSLIDE.R"))
source(paste0(wk.dir, "/lib/plotFeatDensPerRadWindow.R"))
################################################################################
# MAIN CODE * MAIN COhDE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
periphery.lst.len <- length(periphery.lst)
periphery.lst <- rep(periphery.lst, each=length(names(master.lst)))
master.v <- rep(names(master.lst), times=periphery.lst.len)

len <- length(master.v)

for(i in 1:len){
  
  m <- master.v[i]
  lst <- master.lst[[m]]
  periphery <- periphery.lst[[i]]
  
  plotFeatDensPerRadWindow(
    out.dir=out.dir,
    out.name=paste0(out.name, "_", m, "_",
                    ifelse(all(c(1,0)%in%periphery), "periBOTH", paste0("peri", periphery))
                    ),
    DOMXYZRFEATfile=lst$DOMXYZRFEATfile,
    periphery= periphery, 
    foi.list=lst$foi.list,
    feat.type=lst$feat.type,
    foiRData=lst$foiRData,
    # Number of features
    nCPU=nCPU,
    dividingR="SLIDE",
    # Set radius of bin/window
    dr=dr,
    multiplier=lst$multiplier,
    # Cannot generate plot for raw counts because it exceeds the limit
    count.v=lst$count.v,
    # Pick hist.breaks that will not generate 0 densities
    # Should cover all the rVal and may vary between Chrom3D structure
    hist.breaks=hist.breaks,
    regenerateCountData=regenerateCountData,
    regeneratePlotData=regeneratePlotData
  )
  
}
  
# rm(list=ls())

#out.dir=out.dir;
#out.name=paste0(out.name, "_", m);
#DOMXYZRFEATfile=lst$DOMXYZRFEATfile;
#foi.list=lst$foi.list;
#feat.type=lst$feat.type;
#foiRData=lst$foiRData;
## Number of features
#nCPU=nCPU;
## BIN/SLIDE/NONE
#dividingR=dividingR;
## Set radius of bin/window
#dr=dr;
#multiplier=lst$multiplier;
## Cannot generate plot for raw counts because it exceeds the limit
#count.v=lst$count.v;
## Pick hist.breaks that will not generate 0 densities
## Should cover all the rVal and may vary between Chrom3D structure
#hist.breaks=hist.breaks;
#regenerateCountData=regenerateCountData;
#regeneratePlotData=regeneratePlotData