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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/14_ArcPlot"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    #lib = "/t1-data/user/ltamon/DPhil/lib"
    #wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/17_Structure"
    #data.dir= "/t1-data/user/ltamon/Database"
    lib = "/stopgap/sahakyanlab/ltamon/DPhil/lib"
    wk.dir = "/stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/14_ArcPlot"
    data.dir= "/stopgap/sahakyanlab/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_gapVsCp")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(1:22, "X"))
gcb = "min2Mb"
bin.len = 40000
gap.type = "gap.perc" # "gap.bp" | gap.perc
Cp.v = 1:21

plot.type = "box" # hexbin | box | both
# For hexbin plot
cuts=6
bins=60

plotOnly = TRUE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(ggplot2)
library(Hmisc)
library(hexbin)
library(viridis)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/makeHexbinggplot.R"))
### FUNCTION ###################################################################
makebp <- function(out.name=paste0(out.dir, "/", out.name, "_AllLR"),
                   GAP=GAP, xlab=levels(GAP$Cp), ylab=gap.type
){
  
  xlab.len <- length(xlab)
  
  if(gap.type=="gap.perc"){
    max.y <- 100
  } else if(gap.type=="gap.bp"){
    max.y <- ceiling(max(GAP$gap))
  }

  png(file=paste0(out.name, ".png"), width=10, height=10, units="in", res=1200)
  
  # range=0, coef=0 causes the whiskers to extend to the data extremes (so no  
  # data will be outliers)
  bp <- boxplot(gap~Cp, outline=TRUE, data=GAP,  xlab=expression(bold("c"["p"])), 
                ylab=ylab, width=rep(0.4, times=xlab.len), cex.axis=1, col="#FDC776", 
                xaxt="n", range=1.5, ylim=c(0,max.y))
  axis(side=1, at=1:xlab.len, labels=xlab, cex.axis=1)
  mtext(side=3, line=2, cex=0.7, 
        text=paste0(out.name, 
                    "\n whiskers 1.5IQR; Cp=0 includes all LR contacts
                    gap.bp unit is Mb; gap.perc is % of chr length"))
  dev.off()
  
  return(bp$stats)
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
Cp.v <- as.character(sort(unique(Cp.v), decreasing=FALSE))

# Chromosome length file
chrLen.df <- read.table(file=chrLenfile, stringsAsFactors=FALSE, as.is=FALSE, header=TRUE,
                        colClasses=c("character", "integer", "integer"))

GAP <- list()
p.lst <- list()
for(chr in chr.v){
  
  if(plotOnly==FALSE){
    
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    # Gap distance between contacting loci in bp
    #gap.bp <- ((PERSIST.MX$hits$j-PERSIST.MX$hits$i)*bin.len)-bin.len
    gap.bp <- (PERSIST.MX$hits$j-PERSIST.MX$hits$i-1)*bin.len
    gap.perc <- 100*gap.bp/( chrLen.df[chrLen.df$chromosome==chr,"length.bp"] )
    df <- cbind.data.frame(chr=chr, gap.bp=gap.bp, gap.perc=gap.perc, Cp=PERSIST.MX$ntis)
    rm(PERSIST.MX, gap.perc, gap.bp)
    gc()
    
    df$chr <- as.character(df$chr)
    save(df, file=paste0(out.dir, "/", chr, "_", gcb, "_gapVsCp.RData"))
    
  } else {
    load(file=paste0(out.dir, "/", chr, "_", gcb, "_gapVsCp.RData"))
  }
  
  # Hexbin plot
  if( plot.type%in%c("hexbin", "both") ){
    
    p.lst[[chr]] <- makeHexbinggplot(xvar=df$Cp, 
                                     yvar=(df[[gap.type]]), 
                                     bins=bins, 
                                     cuts=cuts,
                                     breaks.y=c(0, 20, 40, 60, 80, 100),
                                     #limits.y=c(0, 100),
                                     xlab=expression(bold("c"["p"])),
                                     ylab=gap.type,
                                     title=paste0(chr, "_", gcb, "_bins=", bins, "_cuts=", cuts),
                                     col=viridis(cuts)
    )[["hexplot"]]
    
    ggsave(filename=paste0(out.dir, "/", chr, "_", gcb, "_bins", bins, "_cuts", 
                           cuts, "_", gap.type, "_gapVsCp.pdf"),
           units="in", width=10, height=10, plot=p.lst[[chr]])
    
  }
  
  # Boxplot data
  if( plot.type%in%c("box", "both") ){
    GAP[[chr]] <- df
  }
  
  rm(df)

  print(paste0(chr, " done!"), quote=FALSE)
  
}

# Hexbin plot
if( plot.type%in%c("hexbin", "both") ){
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=4, ncol=6)
  ggexport(p.arr, height=40, width=60,
           filename=paste0(out.dir, "/", gcb, "_bins", bins, "_cuts", cuts, 
                           "_", gap.type, "_gapVsCp.pdf"))
  
}

#-------------------Boxplot

if( plot.type%in%c("box", "both") ){
  
  GAP <- do.call("rbind", GAP)
  rownames(GAP) <- NULL
  data.table::setnames(x=GAP, old=gap.type, new="gap", skip_absent=FALSE)
  if(gap.type=="gap.bp"){
    GAP$gap <- GAP$gap/1e6
  }
  
  out.name <- paste0(gcb, "_LR_gapdist_", gap.type)
  
  # Boxplot across Cp
  
  GAP$Cp <- factor(as.character(GAP$Cp), levels=c("0",as.character(Cp.v)))
  bp.stat <- makebp(out.name=paste0(out.dir, "/", out.name, "_acrossCp"), GAP=GAP, 
                    xlab=levels(GAP$Cp), ylab=gap.type)
  bp.stat <- bp.stat[,-1] # Remove stat for factor "0", which are all NAs
  colnames(bp.stat) <- as.character(Cp.v)
  rownames(bp.stat) <- c("Q1minus1.5IQR", "25th", "Median", "75th", "Q3plus1.5IQR")
  
  # Boxplot for all LR contacts
  
  GAP$Cp[GAP$Cp%in%Cp.v] <- "0"
  temp <- makebp(out.name=paste0(out.dir, "/", out.name, "_All"), GAP=GAP, 
                 xlab=levels(GAP$Cp), ylab=gap.type)
  bp.stat <- cbind(AllLR=temp[,1], bp.stat)
  rm(temp, GAP)
  write.csv(x=bp.stat, file=paste0(out.dir, "/", out.name, "_stat.csv"), 
            row.names=TRUE, quote=FALSE)
  
}
  
# rm(list=ls()); gc()