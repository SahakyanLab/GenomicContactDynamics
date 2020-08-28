################################################################################
# HicRepeatVsPersist clustering
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "fam" # "fam" | "subfam" | "subfam6"
data.dir = paste0(wk.dir, "/out_HicRepeatHeatmap/", rep.group)
out.dir = paste0(wk.dir, "/out_HicRepeatClustering/", rep.group)
### OTHER SETTINGS #############################################################
# Age rank, ELMTISSDYN identifier
elm.id = "GiorPubl" 
suffix.v = c("GiorPubl_RAW", "GiorPublbinned_NORM1")
gcb = "min2Mb" 
chr = "chrALL"
silhouette = FALSE
clustering = TRUE
numClusters = 3
SEED = 438
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
library(RColorBrewer)
library(factoextra)
source(paste0(lib, "/checkForInfInMx.R"))
## FUNCTION ####################################################################
# Check for NA/NaN/Inf in matrices because fviz_nbclust cannot handle those
#checkForInfiniteVal <- function(mx = ELMTISSDYN.MX.norm){
#  v <- apply( X=mx, MARGIN=1, FUN=function(ntiscol){sum(is.finite(ntiscol))} )
#  return(which(v==FALSE))
#}

HicRepeatCluster <- function(
  out.dir = paste0(wk.dir, "/out_clustering"),
  MX = ELMTISSDYN.MX.norm,
  suffix = "GiorPublbinned_NORM1",
  gcb = "min2Mb", #c("min2Mb", "min05Mb")
  chr = "chrALL",
  silhouette = TRUE,
  clustering = TRUE,
  numClusters = 2,
  seed = 123
){
  coul0 <- brewer.pal(11, "Spectral")
  set.seed(seed)
  
  unique.ntis <- as.numeric(dimnames(MX)[[2]])
  coul <- colorRampPalette(rev(coul0))(length(unique.ntis))
  
  # Check for NA/NaN/Inf in normalized matrix because fviz_nbclust and clust 
  # cannot handle those
  # There are NaNs in normalized matrix because of 0s in raw matrix
  NaN.ind <- checkForInfInMx(mx=MX)
  
  # remove those repeats
  if( length(NaN.ind )!=0 ){
    MX <- MX[-(NaN.ind ),]
    drop <- paste0(chr, ":", paste(names(NaN.ind), collapse=",") )
    cat("These repeats were removed because of having NaN values.", drop)
    write(drop, file=paste0(out.dir, "/", gcb, "_", suffix, 
                            "_removedRepeats.txt"))
  }
  
  out.id <- paste0(chr, "_", gcb, "_", suff, "_", numClusters, "cl")
  
  if(silhouette==TRUE){
    
    # to identify optimal number of clusters
    pdf(file=paste0(out.dir, "/", out.id, "_optimalNumClust_silhMetho.pdf"), 
        height=8, width=8)
    
    silh <- fviz_nbclust(x=MX, FUNcluster=kmeans, method="silhouette",
                        k.max=10, nboot=100, barfill="navy", barcolor="navy",
                        linecolor="navy") 
    print(silh)
    dev.off()
    numClusters <- which(max(silh$data[,"y"])==silh$data[,"y"])
    
    print("Silhouette analyses is DONE!", quote=FALSE)
    
  } # silhouette end
  
  if(clustering==TRUE){
    
    clust <- kmeans(x=MX, centers=numClusters, 
                    iter.max = 1000, nstart = 100, trace=FALSE)
    pdf(file=paste0(out.dir, "/", out.id, "_CLUSTlineplot.pdf"), 
        height=20, width=15)
    par( mfrow=c(3,2), mar=c(5.1, 7, 5, 2.1), mgp=c(4, 1.3, 0) )
    
    for(cl in 1:numClusters){
      elements.inclust <- names(which(clust$cluster==cl))
      addtext <- c(paste0(">Cluster", cl), paste0(elements.inclust,collapse=";"))
      write(x=addtext, file=paste0(out.dir, "/", out.id, "_HiCRepeatClusters_Subfam.txt"), 
            ncolumns=1, sep="\t", append=TRUE)
      y.range <- range(MX[elements.inclust,])
      
      plot(NA, ylim=y.range, xlim=range(unique.ntis),
           col=coul, cex=1.5, pch=21, lwd=3, # ylim=c(0,30)
           xlab="",
           ylab="Scaled contact % with >=1 repeat pair",
           cex.main=1, cex.lab=2.5, cex.axis=2.5,
           main=paste0(chr, "_", gcb, "_", suffix, "_clust", cl, "_", 
                       length(elements.inclust)) )
      mtext(side=1, text=expression("c"["p"]), line=5, cex=2.5)
      #mtext(side=2, text="Y axes title", line=3)
      #mtext(side=3, text="Diagram title", line=1.5)
      for(i in 1:21){abline(v=i, col="grey", lty="dotted")}
      for(elm in elements.inclust){
        lines(x=unique.ntis, y=MX[elm,], col="grey", lwd=3)
      }
      lines(x=unique.ntis, y=clust$centers[cl,], col="navy", lwd=4)
    }
    dev.off()
    
    print("HicRepeatClustering is DONE!", quote=FALSE)
    
  } # clustering end
  
  rm(.Random.seed, envir=globalenv())
 
}
################################################################################
HicRepeatCluster <- cmpfun(HicRepeatCluster, options=list(suppressUndefined=TRUE))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Load ELMTISSDYN
load(file=paste0(data.dir, "/", chr, "_", gcb, "_ElmTissDyn_",
                 elm.id, ".RData"))

for(suff in suffix.v){
  
  flenme <- paste0(out.dir, "/", chr, "_", gcb, "_", suff, "_", numClusters,  
                   "cl_HiCRepeatClusters_Subfam.txt")
  if( file.exists(flenme) ){ file.remove(flenme) }
  
  HicRepeatCluster(
    out.dir=out.dir,
    MX=ELMTISSDYN[[suff]],
    suffix=suff,
    gcb=gcb, 
    chr=chr,
    silhouette=silhouette,
    clustering=clustering,
    numClusters=numClusters,
    seed=SEED
  )
  
}

# rm(list=ls()); gc()
