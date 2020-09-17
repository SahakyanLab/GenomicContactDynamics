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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "subfam" # "fam" | "subfam" | "subfam6"
elm.dir = paste0(wk.dir, "/out_HicRepeatHeatmap/viridis/", rep.group)
out.dir = paste0(wk.dir, "/out_AverageBinTrend")
### OTHER SETTINGS #############################################################
# Age rank, ELMTISSDYN identifier
elm.id = "GiorPubl" 
gcb = "min2Mb" 
chr = "chrALL"
bin.size = 31
silhouette = TRUE
# If numClusters=NULL, numClusters will come from silhouette result.
numClusters = NULL
clustering = TRUE
SEED = 438
hmaptype.v = "norm" # c("raw", "norm", "fc")
# Raw >= 0.1; arrange from lowest to highest; only c("MIRb", "MIR") > 0.5;
# marked red in the plot
red.v = c("MER5A", "L3", "L1M5", "AluJo", "L2", "AluJb", "AluSx", "MIR3", "AluY", 
          "MIR", "MIRb")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
library(factoextra)
library(cluster)
source(paste0(lib, "/doVenn.R"))
## FUNCTION ####################################################################
HicRepeatCluster <- function(
  out.dir = paste0(wk.dir, "/out_clustering"),
  MX = ELMTISSDYN.MX.norm,
  suffix = "GiorPublbinned_NORM1",
  gcb = "min2Mb", 
  chr = "chrALL",
  bin.size = 2,
  silhouette = TRUE,
  clustering = TRUE,
  numClusters = 2,
  seed = 123
){

  unique.ntis <- sort(unique( as.numeric(dimnames(MX)[[2]]) ))
  
  # Check for NA/NaN/Inf in normalized matrix because fviz_nbclust and clust 
  # cannot handle those. There are NaNs in norm/fc matrix because of 0s in 
  # raw matrix.
  drp.rw <- apply( X=MX, MARGIN=1, FUN=function(rw) any(!is.finite(rw)) )
  
  # Remove those repeats
  if( sum(drp.rw)>0 ){
    MX <- MX[!drp.rw,]
    drop <- paste0(chr, ":", paste(names(drp.rw==TRUE), collapse=";") )
    cat("These repeats were removed because of having NaN values.", drop)
    write(drop, file=paste0(out.dir, "/", chr, "_", gcb, "_", suffix, 
                            "_removedRepeats.txt"), append=FALSE)
  }
  
  MX.len <- nrow(MX)
  s.ind <- seq(from=1, to=MX.len, by=31)
  e.ind <- c(s.ind[-1]-1, MX.len)
  bin.num <- unique(c(length(s.ind), length(e.ind)))
  
  for( m in c("kmeans_euclidean", "pam_euclidean", "pam_manhattan") ){
    
    clm <- strsplit(x=m, split="_")[[1]]
    dm <- clm[2]; clm <- clm[1]
    out.id <- paste0("seed_", seed, "_", chr, "_", gcb, "_", suff, 
                     "_binsize", bin.size, "_", clm, "_", dm)
    
    #-------------------
    pdf(file=paste0(out.dir, "/", out.id, "_optimalNumClust_silhMetho.pdf"), 
        height=60, width=10)
    par(mfrow=c(bin.num, 2), mar=c(5.1, 7, 5, 2.1), mgp=c(4, 1.3, 0))
    
    for(b in 1:bin.num){
      
      diss.mx <- dist(x=MX[s.ind[b]:e.ind[b],], method=dm)
      
      #-------------------Identify optimal number of clusters]
      if(silhouette==TRUE){
        
        clustm <- ifelse(clm=="kmeans", clm, "cluster::pam")
        
        set.seed(seed)
        eval(parse(text=paste0(
          'silh <- fviz_nbclust(x=MX[s.ind[b]:e.ind[b],], diss=diss.mx, FUNcluster=', clustm, 
          ', method="silhouette", k.max=10, barfill="navy", barcolor="navy", linecolor="navy")' 
        )))
        #print(silh)
        
        if(is.null(numClusters)){
          numClusters <- as.numeric(silh$data$clusters[max(silh$data$y)==silh$data$y])
        }
        
        print(paste0("Cluster ",b, " ", m, ": Silhouette analyses is DONE!"), quote=FALSE)
        
      } # Silhouette end
      rm(diss.mx)
      #-------------------Cluster
      if(clustering==TRUE){
        
        #out.id <- paste0(out.id, "_", numClusters, "cl")
        
        set.seed(seed)
        if(clm=="kmeans"){
          clust <- kmeans(x=MX[s.ind[b]:e.ind[b],], centers=numClusters, 
                          iter.max=1000, nstart=100, trace=FALSE)
          clust <- list(cluster=clust$cluster, centers=clust$centers)
        } else {
          clust <- cluster::pam(x=MX[s.ind[b]:e.ind[b],], k=numClusters, 
                                diss=FALSE, metric=dm, cluster.only=FALSE)
          clust <- list(cluster=clust$clustering, centers=clust$medoids)
        }
        nmecl.v <- rownames(clust$centers)
        
        for(cl in 1:numClusters){
          
          elements.inclust <- names(which(clust$cluster==cl))
          addtext <- c(paste0(">Bin_", b, "_Cluster", cl, "_center_", nmecl.v[cl]),  
                       paste0(elements.inclust, collapse=";"))
          flenme <- paste0(out.dir, "/", out.id, "_HiCRepeatClusters.txt")
          write(x=addtext, file=flenme, ncolumns=1, sep="\t", append=TRUE)
          
          y.range <- range(MX[elements.inclust,])
          plot(NA, ylim=y.range, xlim=range(unique.ntis), cex=1.5, pch=21, lwd=3,
               xlab="", ylab=paste0(suffix, " contact fr with >=1 repeat pair"),
               cex.main=0.5, cex.lab=2.5, cex.axis=2, xaxt="n",
               main=paste0(out.id, "_Bin", b, "_Cluster",  cl, "_", length(elements.inclust)))
          axis(side=1, at=unique.ntis, cex.axis=1, cex.lab=2.5)
          mtext(side=1, text=expression("c"["p"]), line=5, cex=2.5)
          for(i in 1:21){abline(v=i, col="grey", lty="dotted")}
          y.v <- list()
          for(elm in elements.inclust){
            col <- ifelse(elm%in%red.v, adjustcolor("darkred", alpha=0.5),
                          adjustcolor("grey", alpha=0.5))
            if(elm%in%red.v){ y.v[[elm]] <- MX[elm,ncol(MX)] }
            lines(x=unique.ntis, y=MX[elm,], col=col, lwd=3)
          }
          text(x=21, y=unlist(y.v), labels=names(y.v), cex=0.7, col="darkred")
          lines(x=unique.ntis, y=clust$centers[cl,], col=adjustcolor("navy", alpha=0.7), lwd=4)
          rm(y.v)
        }
        
        print(paste0("Cluster ",b, " ", m, ": HicRepeatClustering is DONE!"), quote=FALSE)
        
      } # Clustering end
      
    } # bin.num for loop end
    
    dev.off()
    
  }

  rm(.Random.seed, envir=globalenv())
 
}
################################################################################
HicRepeatCluster <- cmpfun(HicRepeatCluster, options=list(suppressUndefined=TRUE))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.dir <- paste0(out.dir, "/", rep.group)
if( !dir.exists(out.dir) ){
  dir.create(out.dir)
}

load(file=paste0(elm.dir, "/", chr, "_", gcb, "_ElmTissDyn_", elm.id, ".RData"))

for(suff in hmaptype.v){
  
  flenme.v <- paste0(out.dir, "/seed_", SEED, "_", chr, "_", gcb, "_", suff, "_binsize",
                     bin.size, "_", c("kmeans_euclidean", "pam_euclidean", "pam_manhattan"),  
                     "_HiCRepeatClusters.txt")
  
  if( any(file.exists(flenme.v)) ){ file.remove(flenme.v) }
  
  HicRepeatCluster(
    out.dir=out.dir,
    MX=ELMTISSDYN[[suff]],
    suffix=suff,
    gcb=gcb, 
    chr=chr,
    bin.size=bin.size,
    silhouette=silhouette,
    clustering=clustering,
    numClusters=numClusters,
    seed=SEED
  )
  
}

# rm(list=ls()); gc()
