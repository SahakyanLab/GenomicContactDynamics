################################################################################
# Compare contact complementarity distribution of shortest most dynamic contacts
# with that of the longest most persistent contacts. Take top 5% for each group.
# Calculate 5% using all contacts from all chromosomes (not 5% per chromosome)
# because we are more interested in getting contacts at the extremes of the 
# length distribution.
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
    home.dir = "/project/sahakyanlab/ltamon" 
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
complkmer.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity/out_constraints_GfreeSingleNorm/merged_final")
complalign.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity/out_constraints/merged_final")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/25_LengthDependence")
out.dir = paste0(wk.dir, "/out_gap_extremes")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chrALL" #"chrALL" #paste0("chr", c(1:22, "X"))
compl.type = c("kmer", "Gfree", "sdDifference", "align")[INDREPLACE]
bin.len = 40000
dyn.Cp = 1:3
per.Cp = 19:21
perc.thresh = 5 
mult.dyn.thresh = 0.01 # Multiply perc.thresh with this to get perc.thresh for dynamic set
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/compareTwoDist.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( compl.type %in% c("Gfree", "sdDifference") ){
  obj.id <- "kmer"; val.id <- compl.type; compl.dir <- complkmer.dir
} else if( compl.type %in% c("kmer", "align") ){
  obj.id <- compl.type; val.id <- "CII"
  eval(parse(text=paste0( "compl.dir <- compl", compl.type, ".dir" )))
} 

load(paste0(compl.dir, "/", chr, "_", obj.id, "_", gcb, ".RData"))
colnames(CII.MX)[ which(colnames(CII.MX) == "C||") ] <- "CII"
rownames(CII.MX) <- NULL
CII.MX <- CII.MX[ !( is.na(CII.MX[,"Cp"]) | is.na(CII.MX[,val.id]) ), ]
CII.MX <- CII.MX[ CII.MX[,"Cp"] %in% c(dyn.Cp, per.Cp), ]
CII.MX <- cbind.data.frame(CII.MX , gap=CII.MX[,"j"] - CII.MX[,"i"] - 1)
CII.MX <- CII.MX[, setdiff(colnames(CII.MX), c("i", "j")) ]

is_dyn <- CII.MX[,"Cp"] %in% dyn.Cp
is_per <- CII.MX[,"Cp"] %in% per.Cp

dyn.len <- sum(is_dyn) 
per.len <- sum(is_per) 

dyn.max.val <- sort( CII.MX[is_dyn,"gap"], decreasing=F )[ceiling(dyn.len * perc.thresh * mult.dyn.thresh / 100)]
per.min.val <- sort( CII.MX[is_per,"gap"], decreasing=T )[ceiling(per.len * perc.thresh / 100)]

is_dyn <- is_dyn & CII.MX[,"gap"] <= dyn.max.val
is_per <- is_per & CII.MX[,"gap"] >= per.min.val
CII.MX <- CII.MX[ is_dyn | is_per, ]

dyn.id <- paste0(min(dyn.Cp), "To", max(dyn.Cp))
per.id <- paste0(min(per.Cp), "To", max(per.Cp))

is_dyn <- CII.MX[,"Cp"] %in% dyn.Cp
is_per <- CII.MX[,"Cp"] %in% per.Cp 

CII.MX$Cp[is_dyn] <- dyn.id
CII.MX$Cp[is_per] <- per.id
CII.MX$Cp <- factor(CII.MX$Cp, levels=c(per.id, dyn.id))

#
cols <- colorRampPalette(brewer.pal(n=11, name="Spectral"))(21)
cols <- c(adjustcolor(cols[[1]], 1), adjustcolor(cols[[21]], 1))
names(cols) <- levels(CII.MX$Cp)

out.name <- paste0(gcb, "_", chr, "_", bin.len, "_ceilingPerc.thresh", perc.thresh, "_mult.dyn.thresh", mult.dyn.thresh, 
                   "_", compl.type, "_persistent", per.id, "_dynamic", dyn.id)

title.id <- paste0(out.name, "_gap=j - i - 1 in bins_boxplot set to show outliers\n",
                   "Persistent:_Nij=", per.len, "_NijPlot=", sum(is_per), "_perc.thresh=", perc.thresh,
                   " Dynamic:_Nij=", dyn.len, "_NijPlot=", sum(is_dyn), "_perc.thresh=", perc.thresh * mult.dyn.thresh
                   )

pdf(file=paste0(out.dir, "/", out.name, "_box.pdf"), width=20, height=20)
par(mfrow=c(2,2))

p.lst <- list()
for( val.plot in c(val.id, "gap") ){
  
  TEST <- compareTwoDist(x=CII.MX[is_per,val.plot], y=CII.MX[is_dyn,val.plot])
  

  # Density plot
  
  p.lst[[val.plot]] <- ggplot(data=CII.MX, aes_string(x=val.plot)) +
    geom_density(aes(fill=Cp, col=Cp), size=1.5, alpha=0.7) + 
    scale_color_manual(values=cols) + 
    scale_fill_manual(values=cols) +
    labs(title=paste0(title.id, "\n", TEST$test.id)) +
    bgr2 + 
    theme(plot.title=element_text(size=5))
  
  # Box plot
  
  boxplot(as.formula(paste0(val.plot, "~Cp")), outline=F, data=CII.MX, boxwex=0.6, xlab="Cp",
          ylab=val.plot, main=paste0(title.id, "_outline=F \n", TEST$test.id), col=cols, cex.main=0.5, cex.axis=1)
  
  boxplot(as.formula(paste0(val.plot, "~Cp")), outline=T, data=CII.MX, boxwex=0.6, xlab="Cp",
          ylab=val.plot, main=paste0(title.id, "_outline=T \n", TEST$test.id), col=cols, cex.main=0.5, cex.axis=1)

}

dev.off()

p.arr <- ggarrange(plotlist=p.lst, nrow=1, ncol=2)
ggexport(p.arr, height=10, width=20, units="in", filename=paste0(out.dir, "/", out.name, "_dens.pdf"))

# rm(list=ls()); gc()