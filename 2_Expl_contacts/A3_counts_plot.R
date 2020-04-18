################################################################################
# Bar plot of combined long-range contacts per Cp. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/1_Count_contacts"
    data.dir = "/Users/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# Count of long-range contacts
data.dir = out.dir = paste0(wk.dir, "/out_count")
### OTHER SETTINGS #############################################################
gcb.v = c("min2Mb", "min05Mb")
chr.v <- paste("chr", c(1:22, "X"), sep="") 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(gcb in gcb.v){
  
  LR <- read.csv(file=paste0(data.dir, "/", gcb, "_HiCLRcontactsCount.csv"),
                 header=TRUE, row.names=c(chr.v, "chrALL"))
  colnames(LR) <- c("chr", 1:21)
  
  df <- stack(LR[rownames(LR)=="chrALL",-1])
  colnames(df) <- c("count", "Cp")
  # Turn to percentages
  df$count <- df$count/sum(df$count)*100

  ggplot(data=df, aes(x=Cp, y=count)) +
    geom_col(fill="#55bde6", colour="#55bde6") +
    geom_text(aes(label=round(count, digits=2)), vjust=-0.5, size=5) +
    labs(title=paste0(gcb, "_count_LRcontacts_combinedcelltiss"),
         x=expression(bold("c"["p"])),
         y=expression(bold("% Long-range contacts")),
         colour=expression(bold("c"["p"])),
         fill=expression(bold("c"["p"]))
         ) +
    bgr2
  
  ggsave(filename=paste0(out.dir, "/", gcb, "_count_LRcontacts_allcelltiss.pdf"),
         unit="in", width=12, height=12)
  
}

# rm(list=ls())
