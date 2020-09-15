################################################################################
# Visualise maps
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
Cs.raw.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.disc.kmer.5.dir = CII.cont.kmer.5.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/polished/11_Constraints/out_group"
SIM.3.2.kmer.5.5.dir = paste0(wk.dir, "/sim_3.2")
SIM.4.2.kmer.5.dir = paste0(wk.dir, "/sim_4.2")
out.dir = paste0(wk.dir, "/out_generate_map")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1"
ct = "FC"
# Metric name should match source directory name.
# 5 in "CII.disc.kmer.5" is the cutoff percentage for categorisation. disc means
# discrete (categorised CII), cont means continuouos (orig CII). 
# <CII/SIM>.<disc/cont>.<kmer/align>.<(0,100)>
#metric.v = c("CII.disc.kmer.5", "CII.cont.kmer.5", 
#             "Cp", "Cs.raw", "Cs.norm", 
#             "SIM.disc.kmer.5", "SIM.cont.kmer.5") 
metric.v = c("SIM.4.2.kmer.5", "Cs.norm") 
format = "symmetric" # "symmetric" | "square"
bins.i = 1000:3200
bins.j = 1000:3200
# If scalebr.v==NULL, no scale bar
scalebr.v = c(xmin=1030, xmax=1060, ymin=1030, ymax=1060)
gap.range = c(50,1900)
# Useful in case not whole chr is to be plotted
out.id = "whole_gap50To1900_inclu1000To3200"
res = 300

dropUnwantedContacts = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(reshape)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/getContactDF.R"))
source(paste0(wk.dir, "/lib/processForMap.R"))
source(paste0(wk.dir, "/lib/makeMatrixMap.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste(gcb, chr, ct, format, out.id, sep="_")
print(paste0(out.name, "..."), quote=FALSE)

p.lst <- list()
for(metric in metric.v){
  
  eval(parse(text=paste0(
    'metric.dir <- ', metric, '.dir'
  )))
    
  df <- getContactDF(metric.dir=metric.dir, metric=metric, gcb=gcb, chr=chr, 
                     ct=ct, bins.i=bins.i, bins.j=bins.j, gap.range=gap.range)
  if(dropUnwantedContacts){
    df <- df[!is.infinite(df$value),]
  } else{
    df[is.infinite(df$value),"value"] <- NA 
  }

  p.lst[[metric]] <- makeMatrixMap(df=df, format=format, check.dup=FALSE, metric=metric,
                                   plot.title=paste0(out.name, "_", metric, "_scale=", 
                                                     unname(scalebr.v["xmax"]-scalebr.v["xmin"])
                                                     ), 
                                   scalebr.v=scalebr.v)
  
  rm(df); gc()
  print(paste0(metric, " done!"), quote=FALSE)
  
}

m.v.len <- length(metric.v)
p.arr <- ggarrange(plotlist=p.lst, nrow=1, ncol=m.v.len,
                   legend=NULL)
ggexport(p.arr, height=5*res, width=5*m.v.len*res, res=res,
         filename=paste0(out.dir, "/", out.name, ".png"))
ggexport(p.arr, height=10, width=10*m.v.len, 
         filename=paste0(out.dir, "/", out.name, ".pdf"))

# rm(list=ls()); gc()
