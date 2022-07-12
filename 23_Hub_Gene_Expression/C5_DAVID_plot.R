################################################################################
# Plot DAVID functional term clustering results on a heatmap to give identity
# to selected hubs.
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
hub.dir = paste0(wk.dir, "/out_hubfile/FC")
david.dir = paste0(wk.dir, "/out_DAVID/funxAnnoClustering")
out.dir = paste0(wk.dir, "/out_DAVID_plot")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
hub.id = "FC_topCP3_gapBin50"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(stringr)
library(ComplexHeatmap)
library(grid)
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
  file.id <- paste0(david.dir, "/", hub, "_ftc.txt")
  if( file.exists(file.id) ){
    x <- readLines(con=file.id)
  } else {
    
    print(paste0(file.id, " does not exist."), quote=F)
    FUNX[[hub]] <- cbind.data.frame(hub=hub, matrix(data=NA, nrow=1, ncol=13), 
                                    nhubgenes=nhubgenes, nclust1genes=NA,
                                    score=NA, stringsAsFactors=F)
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
  # Get maximum top 3 terms based on PValue
  topN <- ifelse(nrow(cl1.df)>3, 3, nrow(cl1.df))
  # The lower, the more significant
  cl1.df <- cl1.df[order(cl1.df$PValue, decreasing=F),][1:topN,]
  
  # Count genes in Annotation Cluster 1
  nclust1genes <- unlist(strsplit(x=cl1.df$Genes, split=",", fixed=T))
  nclust1genes <- stringr::str_trim(string=nclust1genes, side="both")
  nclust1genes <- length(unique(nclust1genes))
  
  FUNX[[hub]] <- cbind.data.frame(hub=rep(hub), cl1.df, nhubgenes=rep(nhubgenes), 
                                  nclust1genes=rep(nclust1genes), score=rep(score), 
                                  stringsAsFactors=F)
  
  rm(x, file.id, cl.ind, cl.splt, cl.nme, cl.score, nclust1genes, cl1.df, cl1.start, cl1.end, topN)
  
  print(paste0(hub, " done!"), quote=F)
  
} # HUB for loop end

FUNX <- lapply(X=FUNX, FUN=function(x){
  colnames(x) <- colnames(FUNX$min2Mb_chr1_FC_topCP3_gapBin50_hub1)
  return(x)
})

FUNX <- do.call("rbind", FUNX)

rownames(FUNX) <- NULL
FUNX$hub <- gsub(x=FUNX$hub, pattern=paste0(gcb, "_"), replacement="", fixed=T)
FUNX$hub <- gsub(x=FUNX$hub, pattern=paste0("_", hub.id), replacement="", fixed=T)
FUNX$PValue <- as.numeric(FUNX$PValue)
FUNX <- FUNX[order(FUNX$Term, decreasing=F),]

#save(FUNX, file=paste0(out.dir, "/", out.id, "data.RData"))

# Heatmap matrix
term.v <- unique(FUNX$Term[!is.na(FUNX$Term)])
hub.v <- unique(FUNX$hub)

FUNX.MX <- matrix(data=NA, nrow=length(term.v), ncol=HUB.len, dimnames=list(term.v, hub.v))

for(hub in hub.v){
  
  print(paste0(hub, "..."), quote=F)
  tmp <- FUNX[FUNX$hub==hub & !is.na(FUNX$Term),]
  if(nrow(tmp)==0){
    
    print(paste0(hub, " skipped."), quote=F)
    next
    
  } else {
    FUNX.MX[tmp$Term, tmp$hub] <- tmp$PValue 
  }
  rm(tmp)
  
}

# No p-value cut-off applied because I want to give to identity to all hubs
# even with not significant terms, just note this on the plot. Also the EASE
# score is indicated on the plot, a lower value usually means <0.05 p-value of terms
FUNX.MX[is.na(FUNX.MX)] <- 0
FUNX.MX[FUNX.MX>0 & !is.na(FUNX.MX)] <- 1

if( any(!unique(as.vector(FUNX.MX))%in%c(1, 0)) ){
  stop("Error in heatmap values.")
}

# Add EASE score to hub rownames
tmp <- FUNX[,c("hub", "score")]
tmp <- tmp[!duplicated(tmp),]
tmp <- tmp$score[match(x=colnames(FUNX.MX), table=tmp$hub)]
tmp <- stringr::str_trim(string=format(x=tmp, digits=2, scientific=F), side="both")
colnames(FUNX.MX) <- paste(tmp, colnames(FUNX.MX), sep=" - ")
rm(tmp, FUNX)

h <- ComplexHeatmap::Heatmap(matrix=FUNX.MX, col=c("gray80", "darkred"),
                             cluster_columns=T, cluster_rows=F, 
                             clustering_distance_columns="binary",
                             column_dend_height=unit(50, "mm"),
                             column_names_max_height=unit(15, "cm"),
                             column_names_gp=grid::gpar(fontsize=20),
                             column_title=paste0(out.id, "_DAVID_annotationCluster1_maxTop3TermsHighestPvalue\nnopvaluecutoffapplied_alltermsinred_rownameshasEASEofCluster1"), 
                             column_title_gp=gpar(fontsize=10),
                             
                             row_split=1:nrow(FUNX.MX), 
                             row_names_max_width=unit(30, "cm"),
                             row_names_gp=grid::gpar(fontsize=20))

pdf(file=paste0(out.dir, "/", out.id, "_hmap.pdf"), width=30, height=30)
print(h)
dev.off()

#gplots::heatmap.2(x=FUNX.MX, Rowv=F, Colv=T, dendrogram="column", scale="none",
#                  trace="none", na.rm=T, margins=c(10, 30), col=c("gray80", "darkred"), 
#                  cexRow=1, cexCol=1, rowsep=1:nrow(FUNX.MX), colsep=1:ncol(FUNX.MX),
#                  key=F, distfun=function(x) dist(x, method="binary"),
#                  main=paste0(out.id, "_DAVID_annotationCluster1_\nmaxTop3TermsHighestPvalue")) 

# rm(list=ls()); gc()



