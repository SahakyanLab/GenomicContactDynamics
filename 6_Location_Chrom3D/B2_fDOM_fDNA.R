################################################################################
# Plot the density of the DNA and domains at several radial windows of the model
# using DOMXYZR.DF 
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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/6_Location_Chrom3D"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "H1-hESC_LMNB1_hg38" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
# DOMXYZR.DF directory
data.dir = paste0(wk.dir, "out_AddXYZR")
out.dir = paste0(wk.dir, "/out_fDOM_fDNA/", model.id)
### OTHER SETTINGS #############################################################
ploidy = "haploid"
# For IMR90_LMNB1_GSE49341_hg19 diploid, r can be greater than 6
breaks.v = 0:6
limits.v = c(0,6)
# Widths of radial window to be tested
dr.v = c(0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
plotWhich = c("fbeadDNA") # c("fbeadDNA", "fDOM")
HxbinPlot = FALSE
rwindPlot = TRUE
plotOnly = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/BINorSLIDE.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(model.id, "_", ploidy)
# Load DOMXYZR.DF 
load(paste0(data.dir, "/", out.name, "_domXYZR.RData"))
DOMXYZR.DF$domainlen <- abs( as.numeric(DOMXYZR.DF$end)-as.numeric(DOMXYZR.DF$start) )
rmin <- min(DOMXYZR.DF$radDist, na.rm=TRUE)
rmax <- max(DOMXYZR.DF$radDist, na.rm=TRUE)

#----------HEXBINPLOT----------
if(HxbinPlot==TRUE){
  coul <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  p0 <- ggplot(data=DOMXYZR.DF, aes(x=radDist, y=domainlen) ) +
    geom_hex(bins=70) +
    xlab(label="Radial Position") +
    ylab(label="Domain Length") +
    scale_fill_gradientn(colours = coul(11)) +
    ggtitle(label=out.name) +
    bgr1
  ggsave(file=paste0(out.dir, "/", out.name, "_", 
                     "RadDistVsDomainLen_hexbinplot.jpeg"))
} else {
  print("Not making a hexbin plot.")
}

p.lst <- list()

for(dr in dr.v){
  
  #----------RWINDPLOT-----------
  if(rwindPlot==TRUE){
    
    if(plotOnly==FALSE){
      
      rVal <- BINorSLIDE(numVec=DOMXYZR.DF$radDist, action="SLIDE", 
                         numPoints=1000, dr=dr)
      
      # total number of domains
      if( length(unique(DOMXYZR.DF$id))!=nrow(DOMXYZR.DF) ){
        stop("Bead IDs not unique.")
      } else {
        totDOM <- length(DOMXYZR.DF$id)
        totbeadDNA <- sum(DOMXYZR.DF$domainlen)
      }
      
      fDOM.DNA <- sapply(X=rVal, simplify = FALSE, FUN=function(r.inst){
        #defines window
        ind <- which( DOMXYZR.DF$radDist >= (r.inst-dr) & DOMXYZR.DF$radDist <= (r.inst+dr) )
        #number of domains and sum of bead DNA in the r window
        data.frame(fDOM=length(ind)/totDOM, 
                   fbeadDNA=sum(DOMXYZR.DF[ind, "domainlen"])/totbeadDNA,
                   stringsAsFactors=FALSE) 
      })
      
      fDOM.DNA <- cbind.data.frame(rInstant=rVal, do.call("rbind", fDOM.DNA),
                                   stringsAsFactors=FALSE)
      
      save(fDOM.DNA, file=paste0(out.dir, "/", out.name, "_", paste(plotWhich, collapse="_"),
                            "_dr", dr, ".RData"))
      
    } else {
      load(file=paste0(out.dir, "/", out.name, "_", paste(plotWhich, collapse="_"),
                       "_dr", dr, ".RData"))
    } 
    
    fDOM.DNA <- melt(data=fDOM.DNA, id="rInstant")
    fDOM.DNA <- fDOM.DNA[fDOM.DNA$variable%in%plotWhich, ]
    #fDOM.mean <- mean(as.numeric(fDOM.DNA$fDOM), na.rm=TRUE)
    #fDNA.mean <- mean(as.numeric(fDOM.DNA$fbeadDNA), na.rm=TRUE)
    
    label.v <- c("Domain", "DNA")
    col.v <- c("black", "black")
    names(label.v) <- names(col.v) <- c("fDOM", "fbeadDNA")
    
    p.lst[[as.character(dr)]] <- ggplot(data=as.data.frame(fDOM.DNA), aes(x=rInstant, y=value)) +
      geom_point( aes(colour=factor(x=variable, levels=plotWhich)),
                  size=1) +
      #geom_hline( linetype="dashed", colour="darkred", size=0.9,
      #            aes(yintercept=fDOM.mean) ) +
      #geom_hline( linetype="dotdash", colour="darkblue", size=0.7, 
      #            aes(yintercept=fDNA.mean)) +
      #geom_vline( linetype="dashed", colour="black", size=0.7, 
      #            aes(xintercept=rmin) ) + 
      #geom_vline( linetype="dashed", colour="black", size=0.7, 
      #            aes(xintercept=rmax) ) + 
      scale_x_continuous( name=paste0("r±", dr), 
                          breaks=breaks.v,
                          #breaks=0:max(fDOM.DNA$rInstant, na.rm=TRUE)
                          limits=limits.v 
                          ) +
      scale_colour_manual(values=col.v[plotWhich], 
                          labels=label.v[plotWhich]) +
      labs(title=out.name, y=expression(bold(" Genome fraction")),
           colour="") +
      bgr2 
 
    ggsave(filename=paste0(out.dir, "/", out.name, "_", paste(plotWhich, collapse="_"), 
                           "_dr", dr, ".pdf"), 
           units="in", width=10, height=10, plot=p.lst[[as.character(dr)]])
    
  } else {
    print("Not making a plot based on radial position window.")
  }

}

p.arr <- ggarrange(plotlist=p.lst, nrow=4, ncol=4,
                   legend=NULL)
ggexport(p.arr, width=40, height=40,
         filename=paste0(out.dir, "/", out.name, "_", paste(plotWhich, collapse="_"),
                         ".pdf" ))
# rm(list=ls()); gc()




