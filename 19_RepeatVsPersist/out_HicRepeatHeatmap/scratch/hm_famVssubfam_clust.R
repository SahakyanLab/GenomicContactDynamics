################################################################################
#
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
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
famPth = paste0(wk.dir, "/out_HicRepeatHeatmap/plot_GiorPubl372rankrepFamilies.csv")
subfamPth = paste0(wk.dir, "/out_HicRepeatClustering/subfam/chrALL_min2Mb_GiorPublbinned_NORM1_HiCRepeatClusters_Subfam.txt")
cl.nme = c(`Cluster2`="ClusterI", `Cluster1`="ClusterII")
out.dir = paste0(wk.dir, "/out_HicRepeatHeatmap")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
subfam <- readLines(con=subfamPth)

len <- length(subfam)
cl.ind <- seq(from=1,to=len, by=2)
rep.ind <- seq(from=2,to=len, by=2)

sf <- as.list(subfam[rep.ind])
names(sf) <- gsub(x=subfam[cl.ind], pattern=">", replacement="", fixed=TRUE)

sf <- lapply(X=sf, FUN=function(x){
  strsplit(x=x, split=";")[[1]]
})
sf <- stack(sf)
colnames(sf) <- c("repName", "cluster.sf")

fam <- read.csv(file=famPth, header=TRUE)
x <- merge(x=fam, y=sf, by="repName", all=TRUE)
x <- x[order(x$rank, decreasing=FALSE),]

x$cluster.sf <- cl.nme[as.character(x$cluster.sf)]
  
write.csv(x=x, file=paste0(out.dir, "/hm_famVssubfam_clust.csv"))

# rm(list=ls()); gc()
