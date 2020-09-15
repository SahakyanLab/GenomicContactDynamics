################################################################################
# ROC 
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
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
csv.dir = paste0(wk.dir, "/out_compare")
out.dir = paste0(wk.dir, "/out_cutoffPlotPerMetric")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1"
ct = "FC"
metric.v = c(subj="SIM.4.2.kmer.5", ref="Cs.norm")
out.id = "only_gap50To1900_inclu1000To3200_r2" 

# Closed ranges where most contacts belong; to colour points differently
subj.range = c(0,0.0015)
ref.range = c(0,0.4)

confMxMetric.v = c("TPR", "FPR", "PPV")
confMxMetric.v = c("TPR", "TNR", "PPV", "NPV",  "FNR", "FPR", "FDR", "FOR", "PT",
                   "TS", "ACC", "BA", "F1", "MCC", "FM", "BM", "MK")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(viridis)
library(ggpubr)
source(paste0(wk.dir, "/lib/confusionMxMetric.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
id <- paste(paste(names(metric.v), metric.v, sep="_"), collapse="_")
out.name <- paste(gcb, chr, ct, out.id, id, sep="_"); rm(id)
print(paste0(out.name, "..."), quote=FALSE)

COMPIJMX <- read.csv(file=paste0(csv.dir, "/", out.name, ".csv"),
                     header=TRUE, stringsAsFactors=FALSE)

out.range <- (COMPIJMX$c.offsubj<subj.range[1] | COMPIJMX$c.offsubj>subj.range[2]) |
  (COMPIJMX$c.offref<ref.range[1] | COMPIJMX$c.offref>ref.range[2])
group <- out.range
group[group==TRUE] <- "out-range"
group[group==FALSE] <- "in-range"

COMPIJMX <- COMPIJMX[group=="in-range",]

p.lst <- list()
for(metric in confMxMetric.v){
  
  densval <- confMxMetric(CONFMX=COMPIJMX[,c("TP","FP", "TN", "FN"),],
                          metric=metric)
  # Plot, source: https://www.rayshader.com/reference/plot_gg.html
  p.lst[[metric]] <- ggplot(data=cbind(COMPIJMX[,c("c.offsubj", "c.offref")], densval),
                            aes(x=c.offsubj, y=c.offref, z=densval)) + 
    geom_tile(aes(fill=densval)) +
    geom_contour(color="cyan") +
    scale_fill_viridis() + 
    labs(x="cut-off subj", y="cut-off ref", fill=metric,
         title=paste0(out.name, "_", metric)) +
    bgr2
  
  print(paste0(metric, " done!"), quote=FALSE)
  
}

len <- length(confMxMetric.v)
p.arr <- ggarrange(plotlist=p.lst, nrow=5, ncol=5)
ggexport(p.arr, height=50, width=50,
         filename=paste0(out.dir, "/", out.name, ".pdf" ))

# rm(list=ls()); gc()
