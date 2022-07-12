################################################################################
# Get these hub info for each expression dataset:
# a. A csv file of pertinent number for each hub. 
# b. A csv file of unique genes overlapping with all the identified hubs. 
# c. Per hub, make a heatmap of genes vs. tissue showing the level of expression
# of each hub gene. 
# d. Boxplot of fraction of genes with no data in each tissue (n0-nWithData per tissue) 
# per hub. 
# e. Generate csv file combining pertinent columns from *nvalues.csv. 
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
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
hub.dir = paste0(wk.dir, "/out_hubfile/All_LTr")
acceptable.hub.dir = paste0(hub.dir, "/acceptable_n")
out.dir = paste0(wk.dir, "/out_hubinfo/All_LTr")
### OTHER SETTINGS #############################################################
# Protein-coding gene expression data (bulk RNA-seq)
exprData.dir = paste0(wk.dir, "/out_cleanExprData")
exprData.suffix = "cutoff0_LTr_ALL"
src.id.v = c("data1", "data2")
gcb = "min2Mb"
hub.id = "All_topCP3_gapBin50"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(gplots)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(gcb, "_", hub.id)

FRND <- list()
for(src.id in src.id.v){
  
  HUB <- list.files(path=hub.dir, pattern=hub.id, recursive=F)
  HUB <- HUB[grepl(x=HUB, pattern=gcb)]
  if(length(HUB)==0){
    stop( paste0(src.id, ": No hubs.") ) 
  }
  
  expr.df <- read.csv(file=paste0(exprData.dir, "/expr_", src.id, "_", exprData.suffix, 
                                  ".csv"), header=T, stringsAsFactors=F)
  if( any(duplicated(expr.df$Gene.Name)) ){
    stop( paste0(src.id, ": Duplicated gene name in expression data.") )
  }
  tiss.v <- sort(colnames(expr.df[,!colnames(expr.df)%in%c("chr", "Gene.Name")]))
  
  GENES <- list()
  N <- list()
  for(hub in HUB){
    
    chr <- strsplit(x=hub, split="_", fixed=T)[[1]][2]
    hub <- gsub(x=hub, pattern=".csv", replacement="", fixed=T)
    
    genes <- read.csv(file=paste0(hub.dir, "/", hub, ".csv"), header=T, 
                      stringsAsFactors=F)[,"gene"]
    if( all(is.na(genes)) ){
      print(paste0(src.id, ": No genes for ", hub, "."), quote=F)
      next
    }
    genes <- unlist(strsplit(x=genes, split=";"))
    genes <- unique(genes[!is.na(genes)])
    GENES[[hub]] <- cbind(hub=hub, chr=chr, genes=genes)
    
    n0 <- length(genes)
    N[[hub]] <- c(
      chr=chr,
      n0=n0,
      # The final table already contains genes with available data (non-NA) for at
      # least 1 tissue.
      nWithData=sum(genes%in%expr.df$Gene.Name),
      # Per tissue, number of genes in hub with data (non-NA). 
      apply( X=expr.df[expr.df$Gene.Name%in%genes,tiss.v], MARGIN=2, 
             # NA means no data available
             FUN=function(col) sum(!is.na(col)) )
    )
    tmp <- as.numeric(N[[hub]][tiss.v])
    N[[hub]] <- c(N[[hub]], 
                  Nmean=mean(x=tmp, na.rm=T ), 
                  Nmed=median(x=tmp, na.rm=T ))
    rm(tmp)
    
    #-------------------Boxplot data
    
    f.id <- paste0(src.id, "_", hub)
    FRND[[f.id]] <- data.frame(src.id=src.id, 
                               hub=gsub(hub, pattern=paste0(gcb, "_|", hub.id, "_"), replacement=""), 
                               n0=n0, tissue=tiss.v, 
                               # Fraction of genes with no data per tissue
                               frND=as.numeric(N[[hub]][tiss.v]), 
                               stringsAsFactors=F)
    FRND[[f.id]][,"frND"] <- (n0-FRND[[f.id]][,"frND"])/n0
    
    #-------------------Heatmap
    
    hmap <- merge(x=as.data.frame(genes), y=expr.df[,c("Gene.Name", tiss.v)],
                  all.x=T, by.x="genes", by.y="Gene.Name")
    if(nrow(hmap)!=length(genes)){ stop("Missing genes.") }
    rownames(hmap) <- hmap$genes
    hmap <- data.matrix(hmap[,colnames(hmap)!="genes"])
    
    # Categorise gene expression
    mx <- hmap
    mx[hmap<0.5] <- -2
    mx[hmap>=0.5 & hmap<=10] <- -1
    mx[hmap>10 & hmap<=1000] <- 0
    mx[hmap>1000] <- 1
    rm(hmap); gc()
    # Order mx based on rowSum
    mx <- mx[order(rowSums(mx, na.rm=T), decreasing=T),]
    
    p.title <- paste0("nPertiss=", paste(N[[hub]][-(1:2)], collapse=";"))
    p.title <- paste0(hub, "_", src.id, "\n", "n0=", N[[hub]]["n0"], "_nWithData=", 
                      N[[hub]]["nWithData"], "_", p.title)
    pdf(paste0(out.dir, "/", src.id, "/", hub, "_", src.id, "_hmap.pdf"),
        width=30, height=30)
    
    heatmap.2(x=mx, Rowv=F, Colv=F, dendrogram="none", scale="none",
              main=p.title, margins=c(10,10), trace="none", na.color="#4d494c", 
              xlab=NULL, ylab=NULL, breaks=c(-2,-1.5,-1,-0.5,0,0.5, 1), 
              col=c("gray90", "#2171B5", "#2171B5", "#9ECAE1", "#9ECAE1", "#A50F15"), 
              # Note that rows can be lost when text size is set to be big
              # relative to canvas size (WEIRD)
              key=T, key.title=NA, keysize=0.5, cexRow=0.5, 
              distfun=function(x) dist(x, method="binary"))
    
    dev.off()
    
    rm(genes, mx, n0, f.id, p.title); gc()
    
    print(hub, quote=F)
    
  } # HUB for loop end
  
  #-------------------Save N and genes hub info
  
  GENES <- do.call("rbind", GENES)
  N <- do.call("rbind", N)
  N <- cbind(hub=rownames(N), N)
  # Add column indicating if hub has n=1 for at least 1 tissue
  N <- cbind(N, nEq1=apply(X=N[,-1], MARGIN=1, FUN=function(row){
    return(ifelse(any(row%in%c("0", "1")), "Yes", "None"))
  }))
  
  write.csv(N, paste0(out.dir, "/", out.id, "_", src.id, "_nvalues.csv"),
            row.names=F, quote=F)
  fle.nme <- paste0(out.dir, "/", out.id, "_hubgenes.csv")
  if( !file.exists(fle.nme) ){
    write.csv(GENES, file=fle.nme, row.names=F, quote=F)
  }
  
  rm(GENES, HUB, expr.df, tiss.v, N, fle.nme); gc()
  
} # src.id.v for loop end

#-------------------Boxplot

FRND <- do.call("rbind", FRND)
rownames(FRND) <- NULL
FRND$hub <- gsub(x=FRND$hub, pattern="_", replacement="\n", fixed=TRUE)
FRND$hub <- factor(x=as.character(FRND$hub), 
                   levels=sort(unique(as.character(FRND$hub)), decreasing=F))
FRND$src.id <- factor(x=as.character(FRND$src.id), 
                      levels=sort(src.id.v, decreasing=F))
HUB.len <- length(levels(FRND$hub))

n0.title <- FRND[,c("hub", "n0")]
n0.title <- n0.title[!duplicated(n0.title),]
rownames(n0.title) <- n0.title$hub
n0.title <- paste(x=n0.title[levels(FRND$hub),"n0"], collapse=";")

col.v <- c("gray80", "gray50") #c("aquamarine3", "coral")
leg.id <- paste(paste(levels(FRND$src.id), col.v, sep="-"), collapse=",")

pdf(file=paste0(out.dir, "/", out.id, "_frND_bp.pdf"), 
    width=40, height=10)
par(mar = c(6.1, 6.1, 4.1, 4.1))

boxplot(frND~src.id*hub, outline=T, data=FRND, boxwex=0.6, xlab=NULL, ylab=NULL, cex.axis=2.5, col=col.v, 
        xaxt="n", yaxt="n", main=paste0(out.id, "_legend:", leg.id, "\n", n0.title,
                                        "_Fraction of genes with no data")) 
axis(side=1, at=seq(1.5, HUB.len*2, 2), labels=levels(FRND$hub), cex.axis=1.5, 
     mgp=c(8, 2, 0), lwd.ticks=0)
axis(side=2, las=2, cex.axis=3)
abline(v=seq(0.5, HUB.len*2+0.5, by=2), lty=1, lwd=2, col="gray50") 
abline(h=0.5, lty=2, lwd=5, col="coral")
legend(x="topright", legend=levels(FRND$src.id), col=col.v, lty=1, lwd=5, cex=1,
       bty="n")

dev.off()

# rm(list=ls()); gc()

#-------------------Run from top after placing selecting acceptable hubs
# in designated directory.

col.v <- c("hub", "chr", "n0", "nWithData", "Nmean", "Nmed", "nEq1")
# Combine columns from generated *nvalues.csv
d1 <- read.csv(file=paste0(out.dir, "/", gcb, "_", hub.id, "_data1_nvalues.csv"),
               stringsAsFactors=F, header=T)[,col.v]
d2 <- read.csv(file=paste0(out.dir, "/", gcb, "_", hub.id, "_data2_nvalues.csv"),
               stringsAsFactors=F, header=T)[,col.v]
d <- merge(x=d1, y=d2, by="hub", all=T, sort=T, suffixes=c(".Set1", ".Set2"))

if( nrow(d)!=HUB.len ){
  stop("Error in merging.")
}

accHUB <- list.files(path=acceptable.hub.dir, pattern=hub.id, recursive=F)
accHUB <- gsub(x=accHUB, pattern=".csv", replacement="", fixed=T)

acc.TF <- d$hub%in%accHUB
if( sum(acc.TF)==length(accHUB) ){
  d$`expression analysis` <- "N"
  d$`expression analysis`[acc.TF] <- "Y"
} else {
  stop("Acceptable hubs not all in merged table.")
}

write.csv(d, paste0(out.dir, "/", out.id, "_nvalues.csv"), row.names=F, quote=F)

# rm(list=ls()); gc()

