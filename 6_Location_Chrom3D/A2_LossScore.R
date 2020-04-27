################################################################################
# Plot loss scores from the optimization of Chrom3D models
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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    os = "Mac"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "IMR90_LMNB1_GSE49341_hg19" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
lossfile.dir = paste0(wk.dir, "/0_lossLog/", model.id)
out.dir = paste0(wk.dir, "/out_LossScore")
### OTHER SETTINGS #############################################################
ploidy.v = c("haploid", "diploid")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(ploidy in ploidy.v){
  out.name <- paste0(model.id, "_", ploidy, "_LS")
  lossScorefile <- paste0(lossfile.dir, "/", dir(lossfile.dir, pattern=ploidy)) 
  lossScore.df <- read.table(file=lossScorefile, header=FALSE,
                             col.names=c("nth", "lossScore"))
  lossScore.df$nth1e6 <- lossScore.df$nth/1000000
  lossScore.df$lossScore1e5 <- lossScore.df$lossScore/100000
  save(lossScore.df, file=paste0(out.dir, "/", out.name, ".RData"))
  
  ggplot(data=lossScore.df, aes(x=nth1e6, y=lossScore1e5)) +
    geom_line(size=3, colour="darkred") +
    labs(title=out.name, 
         x=bquote(bold("i, "%*%~10^6~"iteration")),
         y=bquote(bold("Loss score, "%*%~10^5))
    ) +
    bgr2
  ggsave(file=paste0(out.dir, "/", out.name, ".pdf"), width=8, height=8)
}

# rm(list=ls())

