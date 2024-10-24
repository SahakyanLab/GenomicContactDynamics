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
    home.dir = paste0("/Users/ltamon")
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
  } else if (whorunsit == "LiezelCluster"){
    home.dir = paste0("/project/sahakyanlab/ltamon")
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")

rep.group = "subfam" # "fam" | "subfam" | "subfam6"
elm.dir = paste0(wk.dir, "/out_HicRepeatHeatmap/subfam_sumrep_atleast2sumrep") 
out.dir = paste0(wk.dir, "/out_HicRepeatClustering/subfam_sumrep_atleast2sumrep")
# Repeats to mark with red line in final plot. Intended for those repeats with 
# high raw fraction.  
#red.v = c("CR1", "TcMar-Tigger", "ERV1", "ERVL", "hAT-Charlie", "ERVL-MaLR", "L2",
#          "Low_complexity", "MIR", "Simple_repeat", "L1", "Alu")r
red.v = "" # c("MIR")
### OTHER SETTINGS #############################################################
# Age rank, ELMTISSDYN identifier
elm.id = "GiorPubl" 
gcb = "min2Mb" 
chr = "chrALL"
silhouette = TRUE
# If numClusters=NULL, numClusters will come from silhouette result.
numClusters = NULL
clustering = TRUE
SEED = 982
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
    drop <- paste0(chr, ":", paste(names(drp.rw)[unname(drp.rw)], collapse=";"))
    cat("These repeats were removed because of having NaN values.", drop)
    write(drop, file=paste0(out.dir, "/", chr, "_", gcb, "_", suffix, 
                            "_removedRepeats.txt"), append=FALSE)
  }
  
  CLUST <- list()
  for( m in c("kmeans_euclidean", "pam_euclidean", "pam_manhattan") ){
    
    
    clm <- strsplit(x=m, split="_")[[1]]
    dm <- clm[2]; clm <- clm[1]
    out.id <- paste0("seed_", seed, "_", chr, "_", gcb, "_", suff, "_", clm, 
                     "_", dm)
    
    diss.mx <- dist(x=MX, method=dm)
    
    #-------------------Identify optimal number of clusters
    if(silhouette==TRUE){
      
      clustm <- ifelse(clm=="kmeans", clm, "cluster::pam")
      
      pdf(file=paste0(out.dir, "/", out.id, "_optimalNumClust_silhMetho.pdf"), 
          height=8, width=8)
      
      set.seed(seed)
      eval(parse(text=paste0(
        'silh <- fviz_nbclust(x=MX, diss=diss.mx, FUNcluster=', clustm, 
        ', method="silhouette", k.max=10, barfill="navy", barcolor="navy", linecolor="navy")' 
      )))
      print(silh)
      
      dev.off()
      
      rm(diss.mx); gc()
      
      if(is.null(numClusters)){
        numClusters <- as.numeric(silh$data$clusters[max(silh$data$y)==silh$data$y])
      }
      
      print(paste0(m, ": Silhouette analyses is DONE!"), quote=FALSE)
      
    } # Silhouette end
    #-------------------Cluster
    if(clustering==TRUE){
      
      out.id <- paste0(out.id, "_", numClusters, "cl")
      
      set.seed(seed)
      if(clm=="kmeans"){
        clust <- kmeans(x=MX, centers=numClusters, iter.max=1000, nstart=100, 
                        trace=FALSE)
        clust <- list(cluster=clust$cluster, centers=clust$centers)
      } else {
        clust <- cluster::pam(x=MX, k=numClusters, diss=FALSE, metric=dm, 
                              cluster.only=FALSE)
        clust <- list(cluster=clust$clustering, centers=clust$medoids)
      }
      nmecl.v <- rownames(clust$centers)
      
      pdf(file=paste0(out.dir, "/", out.id, "_CLUSTlineplot.pdf"), 
          height=20, width=20)
      par(mfrow=c(2,2), mar=c(5.1, 7, 5, 2.1), mgp=c(4, 1.3, 0))
      
      y.range <- range(MX)
      for(cl in 1:numClusters){
        
        elements.inclust <- names(which(clust$cluster==cl))
        CLUST[[paste0(m, "_cl", cl)]] <- elements.inclust
        addtext <- c(paste0(">Cluster", cl, "_center_", nmecl.v[cl]),  
                     paste0(elements.inclust, collapse=";"))
        write(x=addtext, file=paste0(out.dir, "/", out.id, "_HiCRepeatClusters.txt"), 
              ncolumns=1, sep="\t", append=TRUE)
        
        #y.range <- range(MX[elements.inclust,])
        plot(NA, ylim=y.range, xlim=range(unique.ntis), cex=1.5, pch=21, lwd=3,
             xlab="", ylab=paste0(suffix, " contact fr with >=1 repeat pair"),
             cex.main=1, cex.lab=2.5, cex.axis=2, xaxt="n",
             main=paste0(out.id, "_", length(elements.inclust)))
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
        text(x=21, y=unlist(y.v), labels=names(y.v), cex=0.5, col="darkred")
        lines(x=unique.ntis, y=clust$centers[cl,], col=adjustcolor("navy", alpha=0.7), lwd=4)
        rm(y.v)
        
      } 
    
      dev.off()
      
      print(paste0(m, ": HicRepeatClustering is DONE!"), quote=FALSE)
      
    } # Clustering end
    
  } # Clustering method for loop end
  
  out.id <- paste0("seed_", seed, "_", chr, "_", gcb, "_", suff)
  doVenn(vennlist=CLUST, filename=paste0(out.dir, "/", out.id), saveVenndata=TRUE)
  
  rm(.Random.seed, envir=globalenv())
 
}
################################################################################
HicRepeatCluster <- cmpfun(HicRepeatCluster, options=list(suppressUndefined=TRUE))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
#out.dir <- paste0(out.dir, "/", rep.group)
#if( !dir.exists(out.dir) ){
#  dir.create(out.dir)
#}

load(file=paste0(elm.dir, "/", chr, "_", gcb, "_ElmTissDyn_",
                 elm.id, ".RData"))

for( suff in names(ELMTISSDYN) ){
  
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
