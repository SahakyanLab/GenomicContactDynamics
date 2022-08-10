################################################################################
# Generate .RData containing terms of the most enriched DAVID functional annotation 
# cluster (Annotation Cluster 1) per hub. Use results from functional clustering 
# medium stringency but use custom for those with no results from the medium setting. 
# Append topN (argument) terms to hubsum.csv with highest Benjamini p-adjusted value. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/Users/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
dir.id = "All_LTr"
david.dir = paste0(wk.dir, "/out_DAVID/", dir.id, "/funxAnnoClustering")
hub.dir = paste0(wk.dir, "/out_hubfile/", dir.id)
hubsum.dir = paste0(wk.dir, "/out_hubsummary/", dir.id)
out.dir = paste0(wk.dir, "/out_addDAVIDToHubFile/", dir.id)
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
hub.id = "All_topCP3_gapBin50"
topN = 5
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(stringr)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(gcb, "_", hub.id)

HUB <- list.files(hub.dir, recursive=F)
HUB <- HUB[grepl(x=HUB, pattern=gcb, fixed=T)]
HUB <- HUB[grepl(x=HUB, pattern=hub.id, fixed=T)]
HUB <- gsub(x=HUB, pattern=".csv", replacement="", fixed=T)
HUB <- sort(unique(HUB, decreasing=F))
HUB.len <- length(HUB)

FUNX <- list()
for(hub in HUB){
  
  # Count hub genes
  
  genes <- read.csv(file=paste0(hub.dir, "/", hub, ".csv"), header=T, 
                    stringsAsFactors=F)$gene
  genes <- unlist(strsplit(x=genes, split=";"))
  nhubgenes <- length( unique(genes[!is.na(genes) & genes!=""]) )
  rm(genes)  
  
  # Load Annotation cluster file
  
  file.med <- paste0(david.dir, "/", hub, "_ftc_medium.txt")
  file.med.TF <- file.exists(file.med)
  file.cus <- paste0(david.dir, "/", hub, "_ftc_custom.txt")
  file.cus.TF <- file.exists(file.cus)
  
  if(file.med.TF){
    
    x <- readLines(con=file.med)
    stringency <- "medium"
    
  } else if( !file.med.TF & file.cus.TF ) {
    
    x <- readLines(con=file.cus)
    stringency <- "custom"
    
  } else if( !file.med.TF & !file.cus.TF ){
    
    FUNX[[hub]] <- cbind.data.frame(hub=hub, matrix(data=NA, nrow=1, ncol=13), 
                                    nhubgenes=nhubgenes, nclust1genes=NA,
                                    score=NA, stringency=NA, stringsAsFactors=F)
    next
    
  }
  
  # Indices of annotation cluster names 
  cl.ind <- which(grepl(x=x, pattern="Annotation Cluster", fixed=T))
  
  # All hubs have annotation cluster but check in case.
  if( length(cl.ind)==0 ){
    stop(paste0(hub, ": No annotation clusters in file."))
  }
  
  cl.splt <- strsplit(x=x[cl.ind], split="\t", fixed=T)
  cl.nme <- unlist( lapply(X=cl.splt, FUN=function(x)x[1]) )
  cl.score <- unlist( lapply(X=cl.splt, FUN=function(x)x[2]) )
  
  # EASE of each annotation cluster
  cl.score <- as.numeric(gsub(x=cl.score, pattern="Enrichment Score: ",  
                              replacement="", fixed=T))
  
  # Confirm that Annotation Cluster 1 has highest enrichment score (EASE).
  if( max(cl.score)!=cl.score[1] ){
    stop(paste0(hub, ": ",  cl.nme[1], " does not have the highest EASE."))
  }
  
  # First term of Annotation cluster 1
  cl1.start <- cl.ind[1]+1
  
  # Last term of Annotation cluster 1
  if( length(cl.ind)==1 ){
    cl1.end <- length(x)
  } else {
    cl1.end <- cl.ind[2]-1
  }
  
  #  Enrichment score of Annotation Cluster 1
  score <- cl.score[1]
  
  # Annotation Cluster 1
  if((cl1.end-cl1.start)<1){
    stop(paste0(hub, ": Error in getting indices of Annotation Cluster 1 terms."))
  }
  cl1.df <- read.delim(text=x[cl1.start:cl1.end], sep="\t", stringsAsFactors=F)
  cl1.df <- cl1.df[order(cl1.df$PValue, decreasing=F),]
  
  # Count genes in Annotation Cluster 1
  nclust1genes <- unlist(strsplit(x=cl1.df$Genes, split=",", fixed=T))
  nclust1genes <- stringr::str_trim(string=nclust1genes, side="both")
  nclust1genes <- length(unique(nclust1genes))
  
  FUNX[[hub]] <- cbind.data.frame(hub=rep(hub), cl1.df, nhubgenes=rep(nhubgenes), 
                                  nclust1genes=rep(nclust1genes), score=rep(score), 
                                  stringency=stringency, stringsAsFactors=F)
  
  rm(x, cl.ind, cl.splt, cl.nme, cl.score, nclust1genes, cl1.df, cl1.start, cl1.end)
  
  print(paste0(hub, " done!"), quote=F)
  
} # HUB for loop end

FUNX <- lapply(X=FUNX, FUN=function(x){
  colnames(x) <- colnames(FUNX[[1]])
  return(x)
})

save(FUNX, file=paste0(out.dir, "/", out.id, "_DAVIDannocluster1.RData"))

FUNX <- lapply(X=FUNX, FUN=function(hubtrm){
  
  #hubtrm <- hubtrm[hubtrm$Benjamini < 0.05,]

  topNfin <- ifelse(nrow(hubtrm) >= topN, topN, nrow(hubtrm))
  hubtrm <- hubtrm[order(hubtrm$Benjamini, decreasing=F),][1:topNfin,]
  
  return(
    
    data.frame(
      hub=unique(hubtrm$hub),
      MostEnrichedCluster1Terms=paste(hubtrm$Term, collapse=";"),
      MostEnrichedCluster1Benjaminipadj=paste(x=format(hubtrm$Benjamini, digits=4), collapse=";"),
      MostEnrichedCluster1EASE=unique(hubtrm$score),
      DAVIDClusteringStringency=unique(hubtrm$stringency),
      stringsAsFactors=F 
    )
    
  )
  
})

FUNX <- do.call("rbind", FUNX)

# Append to hubsummary
out.name <- paste0(out.id, "_hubsum_withTop", topN, 
                   "termsOfMostEnrichedDAVIDCluster")

hubsum <- read.csv(file=paste0(hubsum.dir, "/", out.id, "_hubsum.csv"))
FUNX <- FUNX[match(x=hubsum$hub, table=FUNX$hub),]

if( identical(hubsum$hub, FUNX$hub) ){
  
  hubsum <- cbind(hubsum, FUNX[,-1])
  write.csv(x=hubsum, file=paste0(out.dir, "/", out.name, ".csv"), row.names=F)
  
} else {
  stop("Not identical!")
}

# rm(list=ls()); gc()
