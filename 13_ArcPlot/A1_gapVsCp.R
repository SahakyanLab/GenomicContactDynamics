################################################################################
# Plot gap distance between contacting loci as % of chromosome length vs. Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/17_Structure"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/17_Structure"
    data.dir= "/t1-data/user/ltamon/Database"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/17_Structure"
    data.dir= "/home/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_gapVsCp")
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(1:22, "X"))
gcb = "min05Mb"
bin.len = 40000
cuts=6
bins=60
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(Hmisc)
library(hexbin)
library(viridis)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/makeHexbinggplot.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Chromosome length file
chrLen.df <- read.table(file=chrLenfile, as.is=FALSE, header=TRUE,
                        colClasses=c("character", "integer", "integer"))
p.lst <- list()
for(chr in chr.v){
  if(plotOnly==FALSE){
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    # Gap distance between contacting loci in bp
    gap.bp <- ((PERSIST.MX$hits$j-PERSIST.MX$hits$i)*bin.len)-bin.len
    gap.perc <- 100*gap.bp/( chrLen.df[chrLen.df$chromosome==chr,"length.bp"] )
    df <- cbind.data.frame(chr=chr, gap.bp=gap.bp, gap.perc=gap.perc, Cp=PERSIST.MX$ntis)
    rm(PERSIST.MX, gap.perc, gap.bp); gc()
    df$chr <- as.character(df$chr)
    save(df, file=paste0(out.dir, "/", chr, "_", gcb, "_gapVsCp.RData"))
  } else {
    load(file=paste0(out.dir, "/", chr, "_", gcb, "_gapVsCp.RData"))
  }
  p.lst[[chr]] <- makeHexbinggplot(xvar=df$Cp, 
                                   yvar=(df$gap.perc), 
                                   bins=bins, 
                                   cuts=cuts,
                                   breaks.y=c(0, 20, 40, 60, 80, 100),
                                   #limits.y=c(0, 100),
                                   xlab=expression(bold("c"["p"])),
                                   ylab=expression(bold("Gap, % of chr length")),
                                   title=paste0(chr, "_", gcb, "_bins=", bins, "_cuts=", cuts),
                                   col=viridis(cuts)
  )[["hexplot"]]
  ggsave(filename=paste0(out.dir, "/", chr, "_", gcb, "_bins", bins, "_cuts", 
                         cuts, "_gapVsCp.pdf"),
         units="in", width=10, height=10, plot=p.lst[[chr]])
  print(paste0(chr, " done!"), quote=FALSE)
}
p.arr <- ggarrange(plotlist=p.lst, nrow=4, ncol=6)
ggexport(p.arr, height=40, width=60,
         filename=paste0(out.dir, "/", gcb, "_bins", bins, "_cuts", 
                         cuts, "_gapVsCp.pdf"))
# rm(list=ls())