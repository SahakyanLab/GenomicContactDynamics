################################################################################
# Visualise Cs and Cp (side by side) of contacts at region of interest using
# PERSIST.MX. 
# bins.x is the vector of bins represeting region of interest. Take all contacts
# formed by bins in bins.x and then filter those contacts by the desired bins.y.
# if bins.y is NULL, all contacts formed by bins.x are taken. 
# Two formats: a. square (x-axis contains bins.x, y-axis contains bins.y)
# b. symmetric (both axes contains bins.x and bins.y)
# Mac, R/3.5.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/13_CsVsCp"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/13_CsVsCp"
    data.dir = "/t1-data/user/ltamon/Database"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/13_CsVsCp"
    data.dir = "/home/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# Matrix containing contacts and Cs/Cp values (PERSIST.MX or MELT.MX)
# PERSIST.MX in this case
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
melt.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
compl.dir = paste0(wk.dir, "/out_constraints")
out.dir = paste0(wk.dir, "/out_visualiseRegion")
### OTHER SETTINGS #############################################################
gcb = "min2Mb" #"" | "min2Mb" | "min05Mb"
chr = "chr1" # "chr12" | "chr12" | "chr1" 
ct = "FC"
# ID of chosen region
region = "chr1_FC" #"KRAS" #"KRAS_lowCslowCp" #"SGIP1_lowCshighCp"
metric.v = c("Cs", "Cp", "CII")
# CII in CII.MX
type = "kmer" #"kmer"
cutoff = 15
# Vector of bins for x-xis
bins.x = NULL #634:636 #1675:1681 #853:859
# Vector of bins for y-axis
# If bins.y=NULL, all contacts formed by bins.x
bins.y = NULL #NULL #3094:3100 #1727:1733 c(1184, 1208, 1210, 1215, 1223, 1228, 1232)
wholeChr = TRUE
scalebr.v = c(xmin=100, xmax=350, ymin=1, ymax=30)
HiC.res = 4e4L
# squarer = half triangle if plotting wholeChr
format = "symmetric" # "symmetric" | "square" | "none"
outtype = "jpeg" # "jpeg" | "pdf"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
library(RColorBrewer)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/multiplot.R"))
source(paste0(wk.dir, "/lib/hmplot.R"))
source(paste0(wk.dir, "/lib/visualiseBinRegions.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(chr, "_", gcb, "_", ct, "_", format, "_", region)

if( wholeChr & is.null(bins.x) & is.null(bins.y) ){
  
  # Load MELT.MX
  load(file=paste0(melt.dir, "/human_", chr, "_allcontacts.RData"))
  
  # Non-contacts for a given cell type will be included because the whole chr
  # should be plotted
  # upper.tri
 
  MX <- rbind(MELT.MX$upper.tri[,c("i", "j")], MELT.MX$upper.tri.nocontact[,c("i", "j")])
  
  if("Cs"%in%metric.v){
    MX <- cbind(MX,
                Cs=c(MELT.MX$upper.tri[,ct], 
                     rep(0, length(MELT.MX$upper.tri.nocontact[,1]))
                     )
                )
  }
  
  rm(MELT.MX); gc()
  
  if("Cp"%in%metric.v){
    # Load PERSIST.MX (original/scaled)
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    
    MX <- cbind(MX,
                Cp=PERSIST.MX$ntis[match(rownames(MX), rownames(PERSIST.MX$hits))]
                )
  }

} else {
  
  # Load PERSIST.MX (original/scaled)
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  
  MX <- PERSIST.MX$hits[ PERSIST.MX$hits[,ct]!=0 , c("i", "j") ]
  
  if("Cs"%in%metric.v){
    MX <- cbind(MX, Cs=PERSIST.MX$hits[ PERSIST.MX$hits[,ct]!=0,ct ])
  }
  
  if("Cp"%in%metric.v){
    MX <- cbind(MX, Cp=PERSIST.MX$ntis[ PERSIST.MX$hits[,ct]!=0 ])
  }
  
}

if("CII"%in%metric.v){
  # Load CII.MX 
  load(file=paste0(compl.dir, "/", chr, "_", type, "_", gcb, "_grouped_cutoff", cutoff, ".RData"))
  MX <- cbind(MX, CII=CII.MX[rownames(MX), "group"])
  out.name <- paste0(out.name, "_", type, "_cutoff", cutoff)
}

if(format=="symmetric"){
  # Rearrange MX in increasing rownames
  ind <- match(as.numeric(rownames(MX)),
               sort(as.numeric(rownames(MX)), decreasing=FALSE))
  MX[ind,] <- MX
  rownames(MX) <- NULL
}

visualiseBinRegions(
  out.dir=out.dir,
  out.name=out.name,
  # MX should be upper triangle for symmetric to work
  MX=MX,
  metric.v=metric.v, 
  bins.x=bins.x,
  bins.y=bins.y,
  format=format,
  scalebr.v=scalebr.v,
  outtype=outtype
  )

# rm(list=ls())
