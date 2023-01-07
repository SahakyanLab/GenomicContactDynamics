################################################################################
# Plots of insertion sites vs. weigted mean Cp or average Cp
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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/18_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "subfam" # "fam" | "subfam"
bin.size = 31
binrep.dir = paste0(wk.dir, "/out_RepeatOverlapPerBinALL/", rep.group)
meanCp.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binWeightedMeanCp")
rankfile = paste0(wk.dir, "/Repeat_rankingbyAge/rep", rep.group, bin.size, ".csv")
out.dir = paste0(wk.dir, "/out_repeatVsmeanCp")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")

wm = "wmeanCp" # "wmeanCp0" | "wmeanCp"

# Sliding window parameters
d = 1 # Radius of sliding window
numSLwind=10 # Number of sliding (overlapping) windows

# Heatmap parameters
cluster = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(viridis)
library(ComplexHeatmap)
source(paste0(lib, "/splitNumericVector.R"))
source(paste0(lib, "/makebp.R"))
source(paste0(wk.dir, "/lib/myheatmap.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(gcb, "_", rep.group, bin.size, "_", wm, "_d", d, "_numSLwind", 
                 numSLwind, "_heatmapcluster", cluster)

if( file.exists(rankfile) ){
  rnk <- read.csv(file=rankfile, header=TRUE, stringsAsFactors=FALSE)
  if(rep.group=="fam"){
    rnk <- cbind.data.frame(repSubfam=rnk$repFamily, repName=rnk$repFamily, 
                            stringsAsFactors=FALSE)
  }
} else {
  stop("rankfile doesn't exist.")
}
rnk.len <- nrow(rnk)

ITEM <- strsplit(x=rnk$repSubfam, split=";", fixed=TRUE)
names(ITEM) <- rnk$repName

HMAP.MX <- matrix(data=NA, nrow=rnk.len, ncol=numSLwind, 
                  dimnames=list(names(ITEM), as.character(1:numSLwind))
)


pdf(file=paste0(out.dir, "/", out.id, "_bp.pdf"), width=20, height=20)
par(mfrow=c(4,4))

for( i in names(ITEM) ){
  
  item <- ITEM[[i]]
  
  X <- list()
  for(chr in chr.v){
    
    # BINREP.MX
    load(file=paste0(binrep.dir,"/", chr, "_BinRep.RData"))
    bins <- as.numeric(BINREP.MX[,"bins"])
    itemf <- intersect(item, dimnames(BINREP.MX)[[2]])
    if( length(itemf)==0 ){
      print(paste0(i, " ", chr, ": Skipped."))
    }
    
    BINREP.MX <- apply(X=cbind(0, BINREP.MX[,itemf]), MARGIN=1, FUN=sum)
    
    # BINWMEANCP.DF
    load(file=paste0(meanCp.dir, "/", gcb, "_", chr, "_weightedMeanCp.RData"))
    if( !identical( as.numeric(BINWMEANCP.DF$bin), bins ) ){
      stop(paste0(chr, ": Bins not in proper order."))
    }
    
    X[[chr]] <- cbind.data.frame(chr=chr, BINWMEANCP.DF, count=BINREP.MX,
                                 stringsAsFactors=FALSE)
    
    print(paste0(chr, " done!"), quote=FALSE)
    
  }
  X <- do.call("rbind.data.frame", c(X, stringsAsFactors=FALSE))
  
  wmeanCp0.TF <- is.finite(X$wmeanCp0)
  wmeanCp.TF <- is.finite(X$wmeanCp)
  
  #for(wm in c("wmeanCp0", "wmeanCp")){
    
    eval(parse(text=paste0( 'x <- X[', wm, '.TF,]' )))
    x <- x[order(x[[wm]], decreasing=FALSE),]
    tot.bin <- nrow(x)
    
    boundsSL <- c( min(x[[wm]]), ceiling(max(x[[wm]])) )
    slide <- splitNumericVector(x=x[[wm]], d=d, action="SLIDE", numSLwind=numSLwind,
                                boundsSL=boundsSL)
    
    slide$midpoints <- NULL
    x <- lapply(X=slide, FUN=function(ind){ x$count[ind] })
    rm(slide)
    x <- stack(x)
    ext <- paste0("boundsSL_", paste(boundsSL, collapse="-"))
    
    # Number of bins per interval 
    binPerInt <- table(x$ind)[levels(x$ind)]
    binPerInt <- paste(paste(names(binPerInt), binPerInt, sep="="), collapse="_")
    binPerInt <- paste0("Totvalidbins=", tot.bin, "\n", binPerInt)
    
    id <- paste0(out.id, "_item", i)
    makebp(df=x, x="ind", y="values", xlab=paste0(wm, "-SLIDE"), ylab=i, 
           addjitter=FALSE, plot.id=paste0(id, "_pointsMEAN_", binPerInt))
    df <- aggregate(x=x$values, by=list(ind=x$ind), FUN=mean)
    HMAP.MX[as.character(i), as.character(df$ind)] <- df$x
    rm(x, df, tot.bin, boundsSL, ext, binPerInt, id); gc()
    
  #}
  
  rm(X, wmeanCp0.TF, wmeanCp.TF); gc()
  
}
  
dev.off()


pdf(file=paste0(out.dir, "/", out.id, "_heatmap.pdf"), width=20, height=20)

hm <- myheatmap(mx=HMAP.MX, colScheme=viridis::viridis(n=10), rep.group=rep.group,
                mx.nme=paste0("raw\n", rep.group, bin.size), cluster=cluster, 
                at.v=seq(0, max(HMAP.MX), by=4))
print(hm)

# norm

HMAP.MX <- (HMAP.MX-rowMeans(HMAP.MX))/apply(X=HMAP.MX, MARGIN=1, FUN=sd)
at.v <- seq(min(HMAP.MX), max(HMAP.MX), by=0.5)

hm <- myheatmap(mx=HMAP.MX, colScheme=viridis::viridis(n=10), rep.group=rep.group,
                mx.nme=paste0("norm\n", rep.group, bin.size), cluster=cluster, 
                at.v=at.v)
print(hm)

dev.off()

# rm(list=ls()); gc()
 