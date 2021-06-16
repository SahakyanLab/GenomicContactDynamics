################################################################################
# Plot distribution (box or density plot) of contact map values and contact
# probability vs genomic distance trendline
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
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
simmap.dir = paste0(wk.dir, "/sim_maps")
Cs.raw.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.disc.kmer.5.dir = CII.cont.kmer.5.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/polished/11_Constraints/out_group"
out.dir = paste0(wk.dir, "/out_boxplot")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr21"
bin.len = 40000

# Select contact maps

#ct.v = "hg19" #"CTREPLACE"
#metric.v = sort(list.files(path=simmap.dir), decreasing=F)[-(1:3)]
ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
              "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC"))[1:3]
metric.v = "Cs.norm"
# Map id format: <cell/tissue>-<metric name>. Metric name should match source 
# directory name, e.g. for metric name Cs.norm directory is Cs.norm.dir. 
# For c||, <CII>.<disc/cont>.<kmer/align>.<(0,100)>, e.g. 5 in "CII.disc.kmer.5" 
# is the cutoff percentage for categorisation. disc means discrete (categorised CII), 
# cont means continuouos (orig CII). 
map.id.v = paste(ct.v, metric.v, sep="-")

# Filter contacts

# If both incl.bin.x and incl.bin.y lists are NULL, use whole chr. 
incl.x = 'incl.bin.x = NULL'
incl.y = 'incl.bin.y = NULL'
mask.x = 'mask.bin.x = list(3038:6232)' #'mask.bin.x = list(3039:6232)'
mask.y = 'mask.bin.y = list(1:3565)'  #'mask.bin.y = list(1:3563)' 
# If vector gap.range is NULL, no filtering. 
gap.v = 'gap.range = c(50, Inf)'

out.id = "whole_maskMidSquare_gap50up_maskx3038To6232y1To3565_xlim0T06"

plotDist = "dens" # "box" | "dens" | NA
xlim.dens = c(0,6)
plotcvd  = F
# See contactprobVsDistance() re its two arguments below
scale.diag.ind = "51" # Diagonal for scaling; depends on gap.range
n.breaks = 30         # Determines number of distance bins to smoothen curve
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(data.table)
library(reshape2)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/contactprobVsDistance.R"))
source(paste0(wk.dir, "/lib/getmapdir.R"))
source(paste0(wk.dir, "/lib/getContactDF.R"))
source(paste0(wk.dir, "/lib/filterContacts.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name0 <- paste(gcb, chr, out.id, sep="_")

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

if( !is.na(plotDist) & plotDist=="dens" ){
  
  pdf(file=paste0(out.dir, "/", gcb, "_", chr, "_", out.id, "_dens.pdf"), 
      width=25, height=25)
  par(mfrow=c(5,5))
  
} 

if(plotcvd){
  CVD <- list()
}

for(map.id in map.id.v){
  
  ct <- strsplit(x=map.id, split="-", fixed=T)[[1]][1]
  metric <- strsplit(x=map.id, split="-", fixed=T)[[1]][2]
  
  out.name <- paste(out.name0, ct, metric, sep="_")
  print(paste0(out.name, "..."), quote=F)
  
  metric.dir <- getmapdir(metric=metric, simmap.dir=simmap.dir)
  
  # Upper triangle contacts only
  df <- getContactDF(metric.dir=metric.dir, metric=metric, 
                     gcb=gcb, chr=chr, ct=ct, gap.range=gap.range, 
                     incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y, 
                     mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y,
                     chrlen.file=chrlen.file, bin.len=bin.len, invalidij.action=NA)
  df$include <- NULL
  posij.len <- length(df$value)
  
  #-------------------Contact probability vs distance
  
  if(plotcvd){
    
    CVD[[map.id]] <- cbind.data.frame(
      map.id=map.id,
      contactprobVsDistance(df=df, bin.len=bin.len, scale.diag.ind=scale.diag.ind, 
                            n.breaks=n.breaks)$smooth,
      stringsAsFactors=F
    )
    
  }
  
  #-------------------
  
  if( !is.na(plotDist) ){
    
    # Store for boxplot; remove NAs
    v <- sort(df$value, decreasing=FALSE, na.last=NA)
    
    # For box and density plots
    
    v.nZero <- v[v>0]
    
    # Use min and max to bound cut-off ranges.
    max.v <- tail(unique(v), n=2)
    min.v <- min(v[v>0], na.rm=TRUE)
    
    validij.len <- length(v)
    nonNA.perc <- format( x=validij.len/posij.len*100, digits=4 )
    non0.perc <- format( x=sum(v>0)/validij.len*100, digits=4 )
    
    if( (sum(v%in%max.v)!=2) & metric%in%c("Cs.norm", "Cs.raw") ){
      print(paste0(metric, ": Cut-off max checkpoint."))
    }
    
    p.title <- paste0(out.name, "\nmin=", min.v, "_max=", max.v[2],
                      "\n", posij.len, chr, "posij_nonNA", nonNA.perc, "%ofchrposij_nonZero", 
                      non0.perc, "%ofnonNAij")
    
  }
  
  rm(df); gc()
  
  #-------------------Boxplot
  
  if( !is.na(plotDist) & plotDist=="box" ){
    
    png(file=paste0(out.dir, "/", gcb, "_", chr, "_", out.id, "_box.png"), 
        res=500, width=6000, height=1500)
    
    # Boxplot excluding 0s
    if( !grepl(x=metric, pattern="CII.") ){
      
      boxplot(x=v.nZero, outline=FALSE, ylab="Value",
              main=paste0(p.title, "_no0sNoOutliers"), cex.main=0.5)
      boxplot(x=v.nZero, outline=TRUE, ylab="Value",
              main=paste0(p.title, "_no0sWithOutliers"), cex.main=0.5)
      
    }
    
    # Boxplot including all values in v
    boxplot(x=v, outline=FALSE, ylab="Value",
            main=paste0(p.title, "_AllNoOutliers"), cex.main=0.5)
    boxplot(x=v, outline=TRUE, ylab="Value",
            main=paste0(p.title, "_AllWithOutliers"), cex.main=0.5)
    
    dev.off()
    
    rm(v.nZero, p.title, max.v, min.v)
  
  }
 
  #-------------------Density plot, non-zero values
  
  if( !is.na(plotDist) & plotDist=="dens" ){
    
    plot(density(x=v.nZero), xlab="Contact value", main=paste0(p.title, "_no0s"),
         cex.main=0.5, xlim=xlim.dens)
    
    rm(v.nZero, p.title, max.v, min.v)
    
  }
  
  print(paste0(map.id, " done!"), quote=F)
  rm(map.id)
  
} # len for loop end

if( !is.na(plotDist) & plotDist=="dens" ){
  dev.off()
}

if(plotcvd){

  CVD <- do.call("rbind", CVD)
  rownames(CVD) <- NULL
  
  p <- ggplot(data=CVD, aes(x=log10(diag.bp), y=log10(val.ave.scaled)) ) +
    geom_line(aes(colour=map.id)) +
    labs(title=paste0(out.name0,"\n AveContactProbScaledToAveofDiag=", 
                      scale.diag.ind, "_n.breaks=", n.breaks),
         x="log10(genomic distance, (i-j)*bin.len)", 
         y="log10(normalised average contact probability per diagonal index)") + 
    bgr2 + 
    theme(axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          legend.text=element_text(size=10))
  ggsave(filename=paste0(out.dir, "/", out.name0, "_cvd.pdf"), height=10, width=10,
         plot=p)
  
}

# rm(list=ls()); gc()
