################################################################################
# Cf vs. Cp boxplots with contact gap filtering generated per tissue but 
# combining data from all chromosomes. Correlation tests and pairwise comparisons
# of Cp distributions performed. Plots contacts with Cp even though they are 
# not present in given tissue but 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}

data.dir = paste0(home.dir, "/Database")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
lib = paste0(home.dir, "/DPhil/lib")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/4_CsVsCp")
out.dir = paste0(wk.dir, "/out_morePlots")

cscp.id = "HiCNormCs" # rawCs |  HiCNormCs
cscp.dir = paste0(wk.dir, "/out_CsVsCpobj/", cscp.id)
print(paste0("CSCP.MX directory: ", cscp.dir))
### OTHER SETTINGS #############################################################
Cp.v = 1:21
gcb = "min2Mb"
chr = "chrALL"
bin.len = 40000
gap.range.bins.closed = c(50,50) # j - i - 1, closed range, NULL to skip filtering
#ct.v <- c("Lu", "Lu") # REMOVE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#library(vioplot)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/compareManyDist.R"))
source(paste0(lib, "/doCorTest.R"))
source(paste0(lib, "/doVarTest.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(persist.dir, "/chr21_Persist_min2Mb.RData"))
ct.v <- sort(setdiff(colnames(PERSIST.MX$hits), c("i", "j")))

out.name <- paste0(gcb, "_", chr, "_", cscp.id, "VsCp_gaprange_", 
                   gap.range.bins.closed[1], "_", gap.range.bins.closed[2])
print(paste0(out.name, "..."), quote=F)

res = 300
png(file=paste0(out.dir, "/", out.name, "_allct_morePlots.png"),
    height=15*res, width=35*res, res=res)
par(mfrow=c(3,7))

VIOL <- list()
VAR <- COR <- PAIR <- list()
ct.v.len <- length(ct.v)
NPERCP <- matrix(data=NA, nrow=length(Cp.v), ncol=ct.v.len, dimnames=list(as.character(Cp.v), ct.v))

for(ct.ind in 1:ct.v.len){
  
  ct <- ct.v[[ct.ind]]
  
  load(paste0(cscp.dir, "/", chr, "_", gcb, "_", bin.len, "bpbin_", ct, "_CsVsCp.RData"))
  
  if( !is.null(gap.range.bins.closed) ){
    CPCS.MX <- CPCS.MX[ CPCS.MX[,"gap.jminusiminus1.bin"] >= gap.range.bins.closed[1] & 
                        CPCS.MX[,"gap.jminusiminus1.bin"] <= gap.range.bins.closed[2], ]
    gap.range.bins.derived <- range(CPCS.MX[,"gap.jminusiminus1.bin"])
  } else {
    gap.range.bins.derived <- NULL
    print(paste0(ct, ": No gap filtering."))
  }
  
  CPCS.MX <- as.data.frame(CPCS.MX)
  CPCS.MX$Cp <- factor(as.character(CPCS.MX$Cp), levels=as.character(Cp.v))
  
  minval.allCp <- format(min(CPCS.MX$Cs), scientific=3, digits=4)
    
  plot.title <- paste0(out.name, "_", ct, "_minimumCsAllCp", minval.allCp,
                       "_Derivedgaprange_", gap.range.bins.derived[1], "_", gap.range.bins.derived[2])
  
  #vioplot(formula=Cs~Cp, data=CPCS.MX, las = 1, plotCentre="line",
  #        xlab="Cp", ylab="Cf", main=paste0(plot.title, "_horizontalLineIsMedian"), col="#FDC776", cex.main=0.8)
  
  #boxplot(formula=Cs~Cp, data=CPCS.MX, outline=T, 
  #        xlab="Cp", ylab="Cf", main=paste0(plot.title, "_outline=F"), col="#FDC776", cex.main=0.8)
  
  # Boxplot
  
  boxplot(formula=Cs~Cp, data=CPCS.MX, outline=T, main=NULL, col="#FDC776", 
          xaxt="n", yaxt="n", xlab=NULL, ylab=NULL)
  axis(side=1, labels=F, at=as.character(Cp.v))
  axis(side=2, las=2, cex.axis=2)

  # Violin
  
  VIOL[[ct]] <- ggplot(data=CPCS.MX, aes(x=Cp, y=Cs)) +
    geom_violin(scale="width", fill="#FDC776", col="black", trim=T, position=position_dodge(0.5)) +
    stat_boxplot(geom="errorbar", width=0, lty="dashed") + 
    stat_summary(fun="median", pch=1, position=position_dodge(0.5)) +
    labs(x=NULL, y=NULL) + 
    bgr1 + 
    theme(axis.text.x=element_blank(), plot.title=element_text(size=7), aspect.ratio=0.8) 
    
  # Significance 
  
  VAR[[ct]] <- tryCatch(
    { doVarTest(xval=CPCS.MX$Cs, grp=CPCS.MX$Cp, out.dir=out.dir, out.name=paste0(out.name, "_", ct)) },
    error = function(e){ NULL }
  )
  
  COR[[ct]] <- tryCatch(
    {
      doCorTest(xval=as.numeric(as.character(CPCS.MX$Cp)), yval=CPCS.MX$Cs, alt="two.sided", 
                exactpval=F, out.dir=out.dir, out.name=paste0(out.name, "_", ct))
    },
    error = function(e){ NULL }
  )
  
  PAIR[[ct]] <- tryCatch(
    {
      compareManyDist(xval=CPCS.MX$Cs, grp=CPCS.MX$Cp, alt="two.sided", out.dir=out.dir, out.name=paste0(out.name, "_", ct))
    },
    error = function(e){ NULL }
  )

  # N per Cp
  NPERCP[,ct] <- table(CPCS.MX$Cp)[as.character(Cp.v)]
  
  print(paste0(ct, " done!"))
  
  rm(CPCS.MX)
  gc()
  
}

dev.off()

write.csv(NPERCP, file=paste0(out.dir, "/" , out.name, "_numberPointsPerCp.csv"), row.names=T)

save(VAR, file=paste0(out.dir, "/" , out.name, "_varbasedtest.RData"))
save(COR, file=paste0(out.dir, "/" , out.name, "_cortest.RData"))
save(PAIR, file=paste0(out.dir, "/" , out.name, "_pairwisedifftest.RData"))

p.arr <- ggarrange(plotlist=VIOL, nrow=3, ncol=7)
ggexport(p.arr, width=35, height=15, filename=paste0(out.dir, "/", out.name, "_violin.pdf" ))

# rm(list=ls()); gc()

