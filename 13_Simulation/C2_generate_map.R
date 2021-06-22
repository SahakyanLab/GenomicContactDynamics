################################################################################
# Visualise contact maps
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
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
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
simmap.dir = paste0(wk.dir, "/sim_maps")
Cs.raw.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.disc.kmer.5.dir = CII.cont.kmer.5.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/polished/11_Constraints/out_group"
out.dir = paste0(wk.dir, "/out_boxplot")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr21"
bin.len = 40000

# Select contact maps

#ct.v = "hg19" #"CTREPLACE"
#metric.v = sort(list.files(path=simmap.dir), decreasing=F)[-(1:3)]
ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
              "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC"))[1:3]
metric.v = "Cs.norm"
# Map id format: <cell/tissue>-<metric name>. Metric name should match source 
# directory name, e.g. for metric name Cs.norm directory is Cs.norm.dir. 
# For c||, <CII>.<disc/cont>.<kmer/align>.<(0,100)>, e.g. 5 in "CII.disc.kmer.5" 
# is the cutoff percentage for categorisation. disc means discrete (categorised CII), 
# cont means continuouos (orig CII). 
map.id.v = paste(ct.v, metric.v, sep="-")

# Filter contacts

# If both incl.bin.x and incl.bin.y lists are NULL, use whole chr. 
incl.x = 'incl.bin.x = NULL'
incl.y = 'incl.bin.y = NULL'
mask.x = 'mask.bin.x = list(3038:6232)' #'mask.bin.x = list(3039:6232)'
mask.y = 'mask.bin.y = list(1:3565)'  #'mask.bin.y = list(1:3563)' 
# If vector gap.range is NULL, no filtering. 
gap.v = 'gap.range = c(50, Inf)'
# Specify what to do with unwanted contacts; NA - set their values to NA, 
# "drop" - remove them from df, "none" - keep them as it is'
invalidij.action = NA # NA | "drop" | "none"

# Useful in case not whole chr is to be plotted
out.id = "whole_maskMidSquare_gap50up_maskx3038To6232y1To3565"

# Mark bins along x- or/and y-axis
mark.x = c(1, 3038, 3565, 6232)
mark.y = c(1, 3038, 3565, 6232)
limits.x = NULL #c(825, 1000)
limits.y = NULL #c(825, 1000)

# If scalebr.v==NULL, no scale bar
scalebr.v = c(xmin=1, xmax=100, ymin=6183, ymax=6232)
res = 500
format = "symmetric" # symmetric | square
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(data.table)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/getmapdir.R"))
source(paste0(wk.dir, "/lib/getContactDF.R"))
source(paste0(wk.dir, "/lib/processForMap.R"))
source(paste0(wk.dir, "/lib/makeMatrixMap.R"))
source(paste0(wk.dir, "/lib/filterContacts.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name0 <- paste(gcb, chr, out.id, sep="_")

print(incl.x, quote=FALSE) 
print(incl.y, quote=FALSE) 
print(mask.x, quote=FALSE) 
print(mask.y, quote=FALSE) 
print(gap.v, quote=FALSE) 

eval(parse(text=incl.x))
eval(parse(text=incl.y))
eval(parse(text=mask.x))
eval(parse(text=mask.y))
eval(parse(text=gap.v))

p.lst <- list()
len <- length(map.id.v)
for(map.id in map.id.v){
  
  ct <- strsplit(x=map.id, split="-", fixed=T)[[1]][1]
  metric <- strsplit(x=map.id, split="-", fixed=T)[[1]][2]
  metric.dir <- getmapdir(metric=metric, simmap.dir=simmap.dir)
  
  out.name <- paste(out.name0, ct, metric, sep="_")
  print(paste0(out.name, "..."), quote=F)
    
  # Upper triangle contacts only
  df <- getContactDF(metric.dir=metric.dir, metric=metric, 
                     gcb=gcb, chr=chr, ct=ct, gap.range=gap.range, 
                     incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y, 
                     mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y,
                     chrlen.file=chrlen.file, bin.len=bin.len, 
                     invalidij.action=invalidij.action)
  df$include <- NULL

  p.lst[[map.id]] <- makeMatrixMap(df=df, format=format, check.dup=FALSE, metric=metric,
                                   plot.title=paste0(out.name, "_", metric, "_scale=", 
                                                     unname(scalebr.v["xmax"]-scalebr.v["xmin"])
                                                     ), scalebr.v=scalebr.v, 
                                   mark.x=mark.x, mark.y=mark.y, limits.x=limits.x, limits.y=limits.y)
  
  print(paste0(metric, " done!"), quote=FALSE)
  rm(df, metric.dir, ct, metric, out.name); gc()
  
}


m.v.len <- length(map.id.v)
p.arr <- ggarrange(plotlist=p.lst, nrow=1, ncol=m.v.len, legend=NULL)
ggexport(p.arr, height=5*res, width=10*m.v.len*res, res=res,
         filename=paste0(out.dir, "/", out.name, ".png"))
#ggexport(p.arr, height=10, width=20*m.v.len, 
#         filename=paste0(out.dir, "/", out.name, ".pdf"))

# rm(list=ls()); gc()
