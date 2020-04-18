################################################################################
# List down the minimum, maximum, mean lengths of TAD per and for all celltiss. 
# Density plot of distance between TAD boundaries (TAD sizes) called by Schmitt
# et al. 2016, per celltiss. TAD counting starts at 1st TAD boundary. Script also
# includes a check confirming that all lengths of borders are equal to Hi-C 
# resolution, 40 kb. 
# Mac, R/3.5.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/14_TADboundary"
    data.dir = "/Users/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# Bed files directory
TADb.dir = paste0(data.dir, "/funx_data/masterpool")
out.dir = paste0(wk.dir, "/out_TADsizeDist")
### OTHER SETTINGS #############################################################
ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
          "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
HiC.res = 40000L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggpubr)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/plotLengthDist.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
p.lst <- list()
min.v <- max.v <- mean.v <- rep(NA, times=length(ct.v))
names(min.v) <- names(max.v) <- names(mean.v) <- ct.v

for(ct in ct.v){
  
  # Read bed file
  bed <- read.table(file=paste0(TADb.dir, "/ct_", ct, "_foi_TADboundary_desc_Schmitt"),
                    col.names=c("chr", "start", "end"), header=FALSE)
  
  # Confirm that length of borders are all equal to HiC.res
  border.len <- bed$end-bed$start
  
  if(unique(border.len)!=HiC.res){
    stop(paste0(ct, ":Checkpoint 1."))
  }
  
  chr.v <- unique(bed$chr)
  
  len.v <- list()
  for(chr in chr.v){
    
    # Check if TAD boundaries are sorted from upstream to downstream
    
    log <- bed$chr==chr
    
    if(is.unsorted(bed[log, "start"]) | is.unsorted(bed[log, "end"]) ){
      stop("TAD boundaries not sored in increasing order.")
    }
    
    len <- sum(log)
    len.v[[chr]] <- (bed$start[log][-1])-(bed$end[log][-len])
     
  } # chr.v for loop end
  
  len.v <- unlist(len.v)
  min.v[ct] <- min(len.v)
  max.v[ct] <- max(len.v)
  mean.v[ct] <- mean(len.v)
  
  p.lst[[ct]] <- plotLengthDist(vec=(len.v/10e6), vline.v=c(0.5, 2), col="deepskyblue3",
                                label.x=bquote(bold("TAD size, "%*%~10^6)),
                                out.name=paste0("chrALL_", ct, "_TADsizeDist_Schmitt"), 
                                out.dir=out.dir)
  
  print(paste0(ct, " done!"), quote=FALSE)
  
} # ct.v for loop end

# Minimum and maximum lengths of TAD per and for all celltiss 

x <- rbind(All=c(min(min.v), max(max.v), mean(mean.v)),
           cbind(min=min.v, max=max.v, mean=mean.v)
           )

write.table(x, file=paste0(out.dir, "/Schmitt2016_stat_TAD"), col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

p.arr <- ggarrange(plotlist=p.lst, nrow=3, ncol=7,
                   legend=NULL)
ggexport(p.arr, height=30, width=70, 
         filename=paste0(out.dir, "/TADsizeDist_Schmitt.pdf" ))

# rm(list=ls())
