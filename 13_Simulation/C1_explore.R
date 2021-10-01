################################################################################
# Plot distribution (box or density plot) of contact map values and contact
# probability vs genomic distance trendline
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
    CII.dir = "/Users/ltamon/DPhil/GCD_polished/11_Complementarity/out_group"
  } else if(whorunsit == "LiezelCluster"){
    #lib = "/t1-data/user/ltamon/DPhil/lib"
    #data.dir = "/t1-data/user/ltamon/Database"
    #wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
    #CII.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/polished/11_Constraints/out_group"
    lib = "/stopgap/sahakyanlab/ltamon/DPhil/lib"
    data.dir = "/stopgap/sahakyanlab/ltamon/Database"
    wk.dir = "/stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
    CII.dir = "/stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/pending/11_Constraints/out_group"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    data.dir = "/home/ltamon/Database"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
simmap.dir = paste0(wk.dir, "/sim_maps")
Cs.raw.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.cont.kmer.5.dir = CII.cont.align.5.dir = CII.cont.Gfree.5.dir = CII.dir
CII.disc.kmer.5.dir = CII.disc.align.5.dir = CII.dir
out.dir = paste0(wk.dir, "/out_explore/CII_boxplots")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
bin.len = 40000

#-------------------SELECT CONTACT MAPS

chr.v = paste0("chr", c(1:22, "X"))

# Map id format: <cell/tissue>-<metric name>. Metric name should match source 
# directory name, e.g. for metric name Cs.norm directory is Cs.norm.dir. 
# For c||, <CII>.<disc/cont>.<kmer/align>.<(0,100)>, e.g. 5 in "CII.disc.kmer.5" 
# is the cutoff percentage for categorisation. disc means discrete (categorised CII), 
# cont means continuouos (orig CII).

#ct.v = c(rep(x="hg19", times=10),
#         sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
#                "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")))
#metric.v = c("SIM.int.set1.1.5", "SIM.int.set1.2.0", 
#             "SIM.int.set2.1.5", "SIM.int.set2.2.0",
#             "SIM.int.set3.2.10.5.1.5", "SIM.int.set3.2.10.5.2.0",
#             "SIM.non.set3.2.10.5.1.5", "SIM.non.set3.2.10.5.2.0",
#             "SIM.int.set3.2.1.5.old", "SIM.int.set4.2.1.5.old", 
#             rep(x="Cs.norm", times=21))
ct.v = "hg19"
metric.v = "CII.cont.kmer.5"

if( length(ct.v)!=length(metric.v) ){
  stop("Each element of ct.v and metric.v should correspond such that the two 
       vectors have the same length.")
} else {
  
  map.id.v <- expand.grid(chr.v, paste0(ct.v, "-", metric.v))
  map.id.v <- apply(X=map.id.v, MARGIN=1, FUN=paste, collapse="-")
  
}

#-------------------FILTER CONTACTS

# If both incl.bin.x and incl.bin.y lists are NULL, use whole chr. 
# For upper triangle contacts, i --> y and j --> x
incl.x = 'incl.bin.x = NULL'
incl.y = 'incl.bin.y = NULL'
mask.x = 'mask.bin.x = NULL' #list(3038:6232)' #j #'mask.bin.x = list(3039:6232)'
mask.y = 'mask.bin.y = NULL' #list(1:3565)'    #i #'mask.bin.y = list(1:3563)' 
# If vector gap.range is NULL, no filtering. Treated as closed range. 
#gap.v = 'gap.range = c(50, Inf)'
gap.v = 'gap.range = c(50, Inf)'

CsScaling = "none" # "none" | "SD" | "MAD" | "BOXWHISK" | "MEANSD" | "MEDMAD"

#-------------------SET PLOT PARAMETERS

# Output specifications
out.id = "whole_50ToInf_CII.cont.5"
out.id <- paste0(out.id, "_CsScaling", CsScaling)

# Density plot
xlim.dens = c(NA,NA)
ylim.dens = c(NA,NA)

# From https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
base.pal = c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
pal <- colorRampPalette(base.pal)(length(map.id.v))
set.seed(123)
pal <- pal[ sample(x=length(pal), size=length(pal)) ]
#pal <- adjustcolor(col=pal, alpha.f=0.3)
names(pal) <- map.id.v

plotDist = "box" # "box" | "dens" | "dens.combined" | "none"
plotcvd  = F

# CVD plot
#lsize.v <- c(rep(x=2, times=10), rep(x=1, times=21))
#ltype.v <- c(rep(x="dashed", times=10), rep(x="solid", times=21))
#names(lsize.v) <- names(ltype.v) <- map.id.v
# See contactprobVsDistance() re its two arguments below
#scale.diag.ind = NA # Diagonal for scaling (usually min(gap.v)+1); If NA, no scaling
# Determines number of distance bins to smoothen curve
#n.breaks = 30

res = 500
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/simulation_lib/contactprobVsDistance.R"))
source(paste0(lib, "/simulation_lib/getmapdir.R"))
source(paste0(lib, "/simulation_lib/getContactDF.R"))
source(paste0(lib, "/simulation_lib/filterContacts.R"))
source(paste0(lib, "/simulation_lib/scaling.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( any(duplicated(map.id.v)) ){
  
  print("Duplicated map.id in map.id.v. Applying unique().")
  map.id.v <- unique(map.id.v)
  
}

out.name0 <- paste(gcb, out.id, sep="_")

print(incl.x, quote=F) 
print(incl.y, quote=F) 
print(mask.x, quote=F) 
print(mask.y, quote=F) 
print(gap.v, quote=F) 

eval(parse(text=incl.x))
eval(parse(text=incl.y))
eval(parse(text=mask.x))
eval(parse(text=mask.y))
eval(parse(text=gap.v))

rm(incl.x, incl.y, mask.x, mask.y, gap.v)

if( !plotDist%in%c("none", "dens", "dens.combined", "box") ){
  stop('Invalid plotDist argument - "none", "dens" or "box" only.')
}

if(plotDist=="dens"){
  
  pdf(file=paste0(out.dir, "/", gcb, "_", out.id, "_xlim", 
                  xlim.dens[1], "_", xlim.dens[2], "_dens.pdf"), 
      width=25, height=25)
  par(mfrow=c(5,5))
  
}

if(plotDist=="dens.combined"){
  
  pdf(file=paste0(out.dir, "/", gcb, "_", out.id, "_xlim", 
                  xlim.dens[1], "_", xlim.dens[2], "_denscombined.pdf"), 
      width=10, height=10)
  
}

if(plotcvd){
  CVD <- list(bydiag=list(), smooth=list())
}

for(map.id in map.id.v){
  
  chr <- strsplit(x=map.id, split="-", fixed=T)[[1]][1]
  ct <- strsplit(x=map.id, split="-", fixed=T)[[1]][2]
  metric <- strsplit(x=map.id, split="-", fixed=T)[[1]][3]
  
  out.name <- paste(out.name0, chr, ct, metric, sep="_")
  print(paste0(out.name, "..."), quote=F)
  
  metric.dir <- getmapdir(metric=metric, simmap.dir=simmap.dir)
  
  # Upper triangle contacts only
  df <- getContactDF(metric.dir=metric.dir, metric=metric, 
                     gcb=gcb, chr=chr, ct=ct, gap.range=gap.range, 
                     incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y, 
                     mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y,
                     chrlen.file=chrlen.file, bin.len=bin.len, invalidij.action=NA)
  df <- df[,c("i", "j", "value")]
  posij.len <- length(df$value)
  
  if( !grepl(x=metric, pattern="CII.") ){
    plot.TF <- is.finite(df$value) & df$value>0
  } else {
    plot.TF <- is.finite(df$value)
  }
  
  # Scale Cs values if needed
  
  scaleF.id <- NULL
  if( metric%in%c("Cs.raw", "Cs.norm") & CsScaling!="none" ){

    scaled <- scaling(x=df$value[plot.TF], approach=CsScaling)
    df$value[plot.TF] <- scaled$scaled
    
    scaleF.id <- paste(format(x=scaled$scalingFactor, digits=4, scientific=F), collapse=",")
    rm(scaled); gc()
    
  }
  
  #-------------------Contact probability vs distance
  
  if(plotcvd){
  
    CVDobj <- contactprobVsDistance(df=df, bin.len=bin.len, scale.diag.ind=scale.diag.ind, 
                                    smooth.method="BIN", n.breaks=n.breaks)
    CVD$bydiag[[map.id]] <- cbind.data.frame(map.id=map.id, CVDobj$bydiagonal,
                                             stringsAsFactors=F)
    CVD$smooth[[map.id]] <- cbind.data.frame(map.id=map.id, CVDobj$smooth,
                                             stringsAsFactors=F)
    
  }
  
  #-------------------
  
  if( plotDist%in%c("box", "dens", "dens.combined") ){
    
    # Store for boxplot; remove NAs
    v <- sort(df$value, decreasing=F, na.last=NA)
    
    # For box and density plots
    v.nZero <- df$value[plot.TF]
    
    validij.len <- length(v)
    nonNA.perc <- format( x=validij.len/posij.len*100, digits=4 )
    non0.perc <- format( x=length(v.nZero)/validij.len*100, digits=4 )
    
    p.title <- paste0(out.name, "\nscalingFactor", scaleF.id, "_min=", min(v.nZero), "_max=", max(v.nZero),
                      "\n", posij.len, chr, "posij_nonNA", nonNA.perc, "%ofchrposij_nonZero", non0.perc, "%ofnonNAij")
    
  }
  
  rm(df); gc()
  
  #-------------------Boxplot
  
  if(plotDist=="box"){
    
    png(file=paste0(out.dir, "/", gcb, "_", out.name, "_box.png"), 
        res=500, width=6000, height=1500)
    par(mfrow=c(1,4))
    
    
    # Boxplot values in v.nZero
    boxplot(x=v.nZero, outline=FALSE, ylab="Value",
            main=paste0(p.title, "_no0sNoOutliers"), cex.main=0.5)
    boxplot(x=v.nZero, outline=TRUE, ylab="Value",
            main=paste0(p.title, "_no0sWithOutliers"), cex.main=0.5)
    
    # Boxplot values in v
    boxplot(x=v, outline=FALSE, ylab="Value",
            main=paste0(p.title, "_AllNoOutliers"), cex.main=0.5)
    boxplot(x=v, outline=TRUE, ylab="Value",
            main=paste0(p.title, "_AllWithOutliers"), cex.main=0.5)
    
    dev.off()
  
  }
 
  #-------------------Density plot, non-zero values
  
  if( plotDist=="dens" ){
    
    plot(density(x=v.nZero), xlab="Contact value", main=paste0(p.title, "_no0s"),
         cex.main=0.5, xlim=xlim.dens, col=pal[map.id], lty=1, lwd=2)
    legend(x="topright", legend=map.id, col=pal[map.id], lty=1, lwd=2, cex=0.8,
           title=NULL, text.font=4)
      
  }
    
  if( plotDist=="dens.combined" ){
    
    if(map.id==map.id.v[1]){
      
      plot(density(x=v.nZero), xlab="Contact value", main=paste0(p.title, "_no0s"),
           cex.main=0.5, xlim=xlim.dens, ylim=ylim.dens, col=pal[map.id], lty=1, lwd=0.5)
      
    } else {
      lines(density(x=v.nZero), cex.main=0.5, xlim=xlim.dens, col=pal[map.id], lty=1, lwd=0.5)
    }
   
  }
  
  print(paste0(map.id, " done!"), quote=F)
  rm(map.id)
  
} # len for loop end


if(plotDist=="dens.combined"){
  
  # Add a legend to density plot
  legend(x="topright", legend=map.id.v, col=pal, lty=1, lwd=2, cex=0.8,
         title=NULL, text.font=4)
  
}

if( grepl(x=plotDist, pattern="dens") ){
  dev.off()
}

if(plotcvd){
  
  p.lst <- sapply(X=c("bydiag", "smooth"), simplify=F, FUN=function(typ){
    
    cvd <- do.call("rbind", CVD[[typ]])
    rownames(cvd) <- NULL
    cvd$map.id <- as.factor(as.character(cvd$map.id))
    
    p <- ggplot(data=cvd, aes(x=log10(diag.bp), y=log10(val.ave.plot)) ) +
      geom_line(aes(colour=map.id, size=map.id, linetype=map.id)) +
      scale_colour_manual(values=pal[levels(cvd$map.id)]) + 
      scale_linetype_manual(values=ltype.v[levels(cvd$map.id)]) + 
      scale_size_manual(values=lsize.v[levels(cvd$map.id)]) +
      labs(title=paste0(out.name0,"\n AveContactProbScaledToAveofDiag=", 
                        scale.diag.ind, "_n.breaks=", n.breaks, "_", typ),
           x="log10(genomic distance, (i-j)*bin.len)", 
           y="log10(normalised average contact probability per diagonal index)") + 
      bgr2 + 
      theme(axis.title.x=element_text(size=10),
            axis.title.y=element_text(size=10),
            legend.title=element_text(size=10), 
            legend.text=element_text(size=10))
    #ggsave(filename=paste0(out.dir, "/", out.name0, "_nbreaks", n.breaks, "_cvd.pdf"), 
    #       height=10, width=10, plot=p)
    return(p)
    
  })
  
  rm(CVD)
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=1, ncol=2, legend=NULL)

  ggexport(p.arr, filename=paste0(out.dir, "/", out.name0, "_nbreaks", n.breaks, 
                                  "scaleDiagInd", scale.diag.ind, "_cvd.png"),
           height=5*res, width=5*2*res)
  
  ggexport(p.arr, filename=paste0(out.dir, "/", out.name0, "_nbreaks", n.breaks, 
                                  "scaleDiagInd", scale.diag.ind, "_cvd.pdf"),
           height=10, width=20)
  
}

# rm(list=ls()); gc()

dat = data.frame(A=c(2,3,0,1), B=c(1,4, 1,0), C=c(4,0,1,1), D=c(2,0,0,4))
rownames(dat) <- c("SKy", "Ing", "Lowl", "embow")
dat.mat = as.matrix(dat)
heatmap(dat.mat, Colv = NA, Rowv = NA)

myCol_AB <- c("orange", "orangered", "red", "firebrick")
myCol_CD <- c("aquamarine", "chartreuse", "green", "green4")

heatmap(dat.mat, Colv = NA, Rowv = NA, 
       col = ifelse(dat.mat[, 1:2], myCol_AB, myCol_CD))

