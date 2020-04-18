################################################################################
# Length distribution of homer motifs
# Mac, R/3.6.1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/15_Motif"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir <- paste0(wk.dir, "/out_motiflength")
### OTHER SETTINGS #############################################################
motif.dir <- "/Users/ltamon/prog/homer/motifs"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
source(paste0(lib, "/plotLengthDist.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
file.v <- list.files(path=motif.dir, include.dirs=FALSE, recursive=FALSE,
                      pattern="motif")
len <- length(file.v)

motif.len.v <- rep(NA, len)
for(m in 1:len){
  motifseq <- readLines(con=paste0(motif.dir, "/", file.v[m]))
  motifseq <- strsplit(x=motifseq[1], split="\t|>")[[1]][2]
  motif.len.v[m] <- nchar(motifseq)
}

plotLengthDist(df=data.frame(variable="", value=motif.len.v),
               vline.v=mean(motif.len.v),
               col.v="#55bde6",
               out.name="homer_motiflength_vertebra",
               out.dir=out.dir,
               label.x="L"
                           
)

# rm(list=ls())
