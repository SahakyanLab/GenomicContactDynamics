################################################################################
# Heatmap of features from all region sets
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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw_ALL")
foifile = NULL #paste0(wk.dir, "/foifile/foifile_heatmap2")
permtsum.dir = paste0(wk.dir, "/out_summary/feat_732")
out.dir = paste0(wk.dir,"/out_heatmap2/feat_732")
### OTHER SETTINGS #############################################################
permtsum.id = "nperm10000_seed662_mxmskfr0"
permtsum.v = c("Cp21", "CptopCP3", "CpAllCs1perc")
eval.f.v = c("numOlapA", "comOlap")
ct.v = c("hg19", "ESC", "FC", "LC")
pval.cutoff = 0.05
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(gplots) 
library(viridis)
library(reshape2)
source(paste0(lib, "/finaliseFOI.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
hrow <- expand.grid(ct.v, permtsum.v)
hrow <- paste0(hrow$Var2, "_", hrow$Var1); rm(ct.v)

out.name <- paste0(permtsum.id, "_", paste(eval.f.v, collapse="_"))

if(plotOnly==FALSE){
  
  # Features
  foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)
  ct.v <- unlist(
    lapply(X=strsplit(x=foi.v, split="ct_|\\_"), FUN=function(x)x[2]) 
  )
  desc.v <- unlist(
    lapply(X=strsplit(x=foi.v, split="desc_|.bed"), FUN=function(x)x[2])
  )
  foi1.v <- unlist(
    lapply(X=strsplit(x=foi.v, split="_foi_|_desc_"), FUN=function(x)x[2])
  )
  foilab.v <- paste(ct.v, foi1.v, desc.v, sep="_")

  ct.v.len <- length(unique(ct.v))
  foi.v.len <- length(foi.v)
  ps.v.len <- length(permtsum.v)
  eval.f.v.len <- length(eval.f.v)
  
  HMAP.MX <- list()
  for(ps in 1:3){
    
    ps.id <- permtsum.v[ps]
    load(paste0(permtsum.dir, "/", permtsum.id, "_", ps.id , "_permtsum.RData"))
    if( any(!foi.v%in%PERMTSUM$foifile) ){
      stop(paste0(ps.id, ": Not all foi in PERMTSUM."))
    }
    # Filter PERMTSUM based on foi.v
    incl.TF <- PERMTSUM$foifile%in%foi.v
    PERMTSUM <- sapply(X=c("pval", "alt"), simplify=FALSE, FUN=function(x){
      return(PERMTSUM[[x]][incl.TF,])
    })
    
    PERMTSUM <- lapply(X=PERMTSUM, FUN=function(mx){
      mx <- as.matrix(mx[foilab.v,eval.f.v])
      if(eval.f.v.len==1){ dimnames(mx)[[2]] <- eval.f.v }
      return(mx)
    })
    
    alt <- PERMTSUM$pval
    alt[ PERMTSUM$alt=="greater" ] <- 1
    alt[ PERMTSUM$alt=="less" ] <- -1
    alt[PERMTSUM$pval>=0.05] <- 0
    rm(PERMTSUM); gc()
    if(  any(!unique(alt)%in%c(-1,0,1)) ){
      stop("Checkpoint 1.")
    }
    
    if(eval.f.v.len>1){
      alt <- apply(X=alt, MARGIN=1, FUN=function(rw){
        rw.uniq <- unique(rw)
        rw.uniq.len <- length(rw.uniq)
        if(rw.uniq.len==1){
          return(rw.uniq)
        } else if(rw.uniq.len>1){
          return(0)
        } else {
          stop("Checkpoint 3.")
        }
      })
    } else {
      alt <- alt[,1]
    }
    
    HMAP.MX[[ps.id]] <- cbind.data.frame(ct=paste0(ps.id, "_", ct.v), 
                                         feature=paste0(foi1.v, "_", desc.v), 
                                         value=unname(alt), stringsAsFactors=FALSE)
    print(paste0(ps.id, " done!"), quote=FALSE)
    rm(alt, ps.id); gc()
    
  }
  
  HMAP.MX <- do.call("rbind.data.frame", HMAP.MX)
  HMAP.MX <- reshape2::acast(data=HMAP.MX, formula=feature~ct)
  HMAP.MX <- HMAP.MX[,hrow]; rm(hrow)
  
  save(HMAP.MX, file=paste0(out.dir, "/", out.name, "_hmap_orig.RData"))
  write.csv(HMAP.MX, file=paste0(out.dir, "/", out.name, "_hmap.csv"),
            row.names=TRUE)
  
} else {
  load(file=paste0(out.dir, "/", out.name, "_hmap.RData"))
}

pdf(paste0(out.dir, "/", out.name, "_hmap.pdf"), width=50, height=50)
HMAP.MX[is.na(HMAP.MX)] <- 0.1
x <- heatmap.2(x=HMAP.MX, Rowv=TRUE, Colv=FALSE,
               distfun=function(x) dist(x, method="euclidean"),
               dendrogram="row", scale="none", na.rm=TRUE,
               breaks=c(-1,-0.5,0,0.1,1), col=c(viridis(n=3)[1], "gray90", "gray70", viridis(n=3)[2]), 
               colsep=1:ncol(HMAP.MX), na.color="gray50", trace="none",
               margins=c(10,10), key=TRUE, key.title=NA, keysize=0.5,
               cexCol=1, cexRow=0.2)
print(x)
dev.off()

write.csv(HMAP.MX[rev(x$rowInd),x$colInd], file=paste0(out.dir, "/", out.name, "_hmap_dendoOrder.csv"),
          row.names=TRUE)


# rm(list=ls()); gc()



