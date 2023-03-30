################################################################################
# Function to generate correlation coefficient vs. Cp dimension
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(RColorBrewer)
# source(paste0(wk.dir, "/lib/compareTwoDist.R"))
### FUNCTION ###################################################################
makeGeneExprVsGroupPlot <- function(CORCP.DF, anno.file, exprData.dir, src.id, 
                                    gene.id, expr.cutoff, out.dir, out.name,
                                    removeOutliers, addmean){
  
  g.ind.lst <- by(data=CORCP.DF[,c("g1.ind", "g2.ind")], 
                  INDICES=list(CORCP.DF$group), FUN=function(df) unique(unlist(df, use.names=F)) ) 

  anno.v <- read.delim(file=anno.file, header=T, stringsAsFactors=F, sep="\t")[,"name2"]
  exprData.file <- paste0(exprData.dir, "/expr_", src.id ,"_cutoff0_", gene.id, ".csv")
  expr.df <- read.csv(file=exprData.file, header=T, stringsAsFactors=F, row.names="Gene.Name")
  expr.df <- expr.df[,colnames(expr.df)!="chr"]
  
  g.genes.lst <- sapply(X=names(g.ind.lst), USE.NAMES=T, simplify=F, FUN=function(grp.nme){
    expr.genes <- intersect( anno.v[ g.ind.lst[[grp.nme]] ], rownames(expr.df) )
  })
  
  ngenesprgrp <- lengths(g.genes.lst)
  ngenesprgrp.id <- paste(paste(names(ngenesprgrp), ngenesprgrp, sep="="), collapse=";")

  rm(g.ind.lst, ngenesprgrp)
  
  g.expr.lst <- sapply(X=names(g.genes.lst), USE.NAMES=T, simplify=F, FUN=function(grp.nme){
    unlist(expr.df[ rownames(expr.df) %in% g.genes.lst[[grp.nme]], ], use.names=F)
  })
  
  #
  GRPEXPR.DF <- stack(g.expr.lst)
  GRPEXPR.DF <- na.omit(GRPEXPR.DF)
  
  GRPEXPR.DF$values[GRPEXPR.DF$values < expr.cutoff] <- 0
  
  if(removeOutliers){
    
    # Remove outliers 
    ind.v <- levels(GRPEXPR.DF$ind)
    box.outl <- matrix(data=NA, nrow=length(ind.v), ncol=7,
                       dimnames=list(ind.v, c("group", "lower.whisk", "upper.whisk",
                                              "val.count", "outl.count", "outl.min", "outl.max")))
    box.outl <- as.data.frame(box.outl)
    
    for(ind in ind.v){
      
      ind.TF <- as.character(GRPEXPR.DF$ind)==ind
      bp.stats <- boxplot.stats(x=GRPEXPR.DF$values[ind.TF], coef=1.5)$stats
      outl.TF <- ind.TF & GRPEXPR.DF$values > bp.stats[5]
      
      box.outl[ind, "group"] <- ind
      box.outl[ind, -1] <- c(bp.stats[1], bp.stats[5],
                             sum(ind.TF), sum(outl.TF),
                             range(GRPEXPR.DF$values[outl.TF])
      )
      
      GRPEXPR.DF <- GRPEXPR.DF[!outl.TF,]
      
    }
    
    if( identical(rownames(box.outl), box.outl$group) ){
      write.csv(box.outl, row.names=F, 
                file=paste0(out.dir, "/", out.name, "_exprval_dynamic-persistentAndNonhub0-Hub1_outliers.csv"))
    }
    
  }
  
  # Compare
  compareManyDist(xval=GRPEXPR.DF$values, grp=GRPEXPR.DF$ind, out.dir=out.dir, 
                  out.name=paste0(out.name, "_forGeneExprValues"))
  
  cols <- colorRampPalette(brewer.pal(n=11, name="Spectral"))(21)
  cols <- c(adjustcolor(cols[[1]], 1), adjustcolor(cols[[1]], 1), 
            adjustcolor(cols[[21]], 1), adjustcolor(cols[[21]], 1))
  
  cols.fill <- cols
  cols.fill[c(1,3)] <- "gray80"
  
  # Boxplot
  
  GRPEXPR.DF$values <- log10(GRPEXPR.DF$values)
  
  png(filename=paste0(out.dir, "/", out.name, "_exprval_dynamic-persistentAndNonhub0-Hub1.png"), 
      height=300*10, width=300*10, res=300)
  
  boxplot(formula=values~ind, data=GRPEXPR.DF, outline=T, xlab=NULL, range=0,
          main=paste0(out.name, "\nNgenes: ", ngenesprgrp.id), lwd=4,
          ylab="expression in TPM; all tissues", #ylim=c(0,60), 
          border=cols, col=adjustcolor(cols.fill, 0.7), cex.main=0.8)
  
  if(addmean){
    df.mean <- stack(by(data=GRPEXPR.DF$values, INDICES=GRPEXPR.DF$ind, FUN=mean, na.rm=F))
    points(x=df.mean$ind, y=df.mean$values, cex=4, pch=15)
    rm(df.mean)
  }
  
  dev.off()

  # Violin plot
  #p <- makeViolinPlot(df=GRPEXPR.DF, x.nme="ind", y.nme="values", sd.mult=1, 
  #                    fill.nme="ind", fill.cols=rep("gray90", 4), fill.legend="none",
  #                    plot.title=out.name, ylim.val=c(0,60), showOutlier=F)
  
  #ggsave(filename=paste0(out.dir, "/", out.name, "_exprval_dynamic-persistentAndNonhub0-Hub1_violin.pdf"),
  #       plot=p, height=10, width=10)
  
}

# rm(list=ls()); gc()