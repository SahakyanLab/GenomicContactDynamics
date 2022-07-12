################################################################################
# Compare different expression-based test statistic (mean, sd, variance, sd/mean) 
# of hubs to bootstrapped estimates by making various plots. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
exprData.dir = paste0(wk.dir, "/out_cleanExprData")
exprData.suffix = "cutoff0_LTr_ALL"
boot.dir = paste0(wk.dir, "/out_nonHubBoot/All_LTr/exprCutoff_05")
hub.dir = paste0(wk.dir, "/out_hubfile/All_LTr/acceptable_n")
out.dir = paste0(wk.dir, "/out_hubVsNonhub/All_LTr/exprCutoff_05_acchubs")
### OTHER SETTINGS #############################################################
src.id = "data2"
gcb = "min2Mb"
hub.id = "All_topCP3_gapBin50"
expr.cutoff = 0.5
# BOOTSTRAP parameters
iter = 10000
seed0 = 763

showOutlier = F
col.v <- c(`-1`="mediumpurple3", `0`="tomato", `1`="#51C56AFF")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
library(gplots)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
exprDataPath <- paste0(exprData.dir, "/expr_", src.id ,"_", exprData.suffix, ".csv")
out.id <- paste0(gcb, "_", hub.id, "_", src.id)
HUB <- list.files(hub.dir, recursive=F)

HUB <- HUB[grepl(x=HUB, pattern=gcb)]
HUB <- HUB[grepl(x=HUB, pattern=hub.id)]
HUB.len <- length(HUB)

expr.df <- read.csv(file=exprDataPath, header=T, stringsAsFactors=F)
tiss.v <- colnames(expr.df[,!colnames(expr.df)%in%c("chr", "Gene.Name")])
tiss.v.len <- length(tiss.v)

hmap <- list()
bp <- list()
p.lst <- list()
for(h in 1:HUB.len){
  
  hub <- strsplit(x=HUB[h], split=".csv", fixed=T)[[1]][1]
  
  #-------------------Hub genes
  
  genes <- read.csv(file=paste0(hub.dir, "/", hub, ".csv"), header=T, 
                    stringsAsFactors=F)[,"gene"]
  if( all(is.na(genes)) ){
    print(paste0("No gene for ", hub, "."), quote=F)
    next
  }
  genes <- unlist(strsplit(x=genes, split=";"))
  genes <- unique(genes[!is.na(genes)])
  n0 <- length(genes)
  nWithData <- sum(genes%in%expr.df$Gene.Name)
  
  #-------------------Hub data
  
  hub.stat <- apply(X=expr.df[expr.df$Gene.Name%in%genes,tiss.v], MARGIN=2, 
                    FUN=function(tiss.val){
                      
    tiss.val <- tiss.val[!is.na(tiss.val)]
    tiss.val[tiss.val<expr.cutoff] <- 0
    n <- length(tiss.val)
    stat <- c(nPerTiss=n, MEAN=mean(tiss.val), MEDIAN=median(tiss.val), SD=sd(tiss.val), VAR=var(tiss.val),
              # Based on https://www.ebi.ac.uk/gxa/help/index.html
              NE=sum(tiss.val==0)/n, 
              LE=sum(tiss.val>=expr.cutoff & tiss.val<=10)/n,
              ME=sum(tiss.val>10 & tiss.val<=1000)/n,
              HE=sum(tiss.val>1000)/n)
    stat <- c( stat, SDdivMEAN=as.numeric(stat["SD"]/stat["MEAN"]) )
    # SDdivMEAN NaN (when MEAN=SD=0), assign to 0 because they
    # should not be excluded since they are data. 0 because
    # when all genes are not expressed then there is no variation.
    stat[is.nan(stat)] <- 0
    # stat can have NAs (e.g. sd and var of 1 value is NA, consequently
    # when sd is NA, SDdivMEAN is NA)
    if( any(is.nan(stat)) ){
      stop(paste(hub, ": NaN hub stat."))
    }
    return(stat)
    
  })
  
  nPerTiss <- hub.stat["nPerTiss",]
  hub.stat <- hub.stat[-1,]
  rm(genes)
 
  #-------------------Bootstrapped data 
  
  load(file=paste0(boot.dir, "/", hub, "_", src.id, "_seed", seed0, "_iter", 
                   iter, ".RData"))
  if( !identical(nPerTiss, BOOT.MX[,"n"]) ){
    stop(paste0(hub, ": n not identical."))
  }
  BOOT.MX <- BOOT.MX[,colnames(BOOT.MX)!="n"]
  # Make sure correspondence with hub.stat
  non.stat <- (t(BOOT.MX))[,tiss.v]
  rm(BOOT.MX)
  se.TF <- grepl(x=rownames(non.stat), pattern="sd.", fixed=T, ignore.case=F)
  # Separate stat values from boot se(standard error) or sd
  non.se <- non.stat[se.TF,]
  non.stat <- non.stat[!se.TF,]
  if( !identical(rownames(non.stat), rownames(hub.stat)) ){
    stop(paste0(hub, ": A. Order of stats different."))
  }  
  if( !identical(rownames(non.se), paste0("sd.", rownames(hub.stat))) ){
    stop(paste0(hub, ": B. Order of stats different."))
  }
  
  #-------------------Heatmap data
  
  sigless <- hub.stat < non.stat-non.se
  siggreat <- hub.stat > non.stat+non.se
  s <- sigless | siggreat
  
  hm <- hub.stat
  hm[!s] <- 0
  hm[s & sigless] <- -1
  hm[s & siggreat] <- 1
  hmap[[hub]] <- data.frame(h.id=NA, stat=rownames(hm), hm, stringsAsFactors=F,
                            row.names=NULL)
  rm(hm, sigless, siggreat, s); gc()
  
  # Heatmap name for hub
  h.id <- gsub(x=hub, pattern=paste0(gcb, "_|", paste0("_", hub.id)), replacement="")
  nAve <- format(mean(nPerTiss), digits=4)
  hmap[[hub]][["h.id"]] <- paste0(h.id, " n0", n0, "_nWD", nWithData, "_nAve", nAve)
  
  #-------------------Individual hub plot
  
  # Combine hub and non-hub stats in a melted dataframe for plotting
  hub.stat <- cbind(type="hub", reshape2::melt(hub.stat), se=0)
  non.stat <- cbind(type="non", reshape2::melt(non.stat), 
                    se=reshape2::melt(non.se)[,"value"])
  df <- rbind(non.stat, hub.stat)
  colnames(df) <- c("type", "stat", "tissue", "value", "se")
  df$type <- factor(df$type, levels=c("non","hub"))
  rm(hub.stat, non.stat, non.se); gc()
  
  # Plot
  nPertiss.title <- paste0("nPertiss=", paste(nPerTiss, collapse=";"))
  
  for( stat in unique(df$stat) ){
    
    p.nme <- paste0(".", stat, "_", hub)
    
    p.lst[[p.nme]] <- ggplot(data=df[df$stat==stat,], aes(x=tissue, y=value)) +
      geom_errorbar(aes(ymin=value-se, ymax=value+se, color=type), width=.2) +
      geom_point(aes(color=type), size=4, na.rm=T) +
      scale_color_manual(values=c("gray60", "#3571cc")) + 
      labs(x=NULL, y=stat,
           title=paste0(strsplit(x=hub, split=".csv")[[1]], "_", stat, "\n",
                        "n0=", n0, "_nWithData=", nWithData, "_", nPertiss.title)
      ) + 
      bgr2 +
      theme(axis.text.x=element_text(size=10, colour="black", face="bold", angle=90,
                                     hjust=1),
            plot.title=element_text(size=5)) 
    
  }
  
  #-------------------Boxplot data
  
  bp[[hub]] <- cbind.data.frame(hub=h.id, nAve=nAve, df, stringsAsFactors=F)
  
  rm(df, h.id); gc()
  print(paste0(hub, " done!"), quote=F)
  
} # HUB.len for loop end

hmap <- do.call("rbind.data.frame", c(hmap, stringsAsFactors=F))

bp <- do.call("rbind.data.frame", c(bp, stringsAsFactors=F))
bp$type <- as.character(bp$type)
bp$stat <- as.character(bp$stat)
bp$tissue <- as.character(bp$tissue)
#bp$type <- factor(as.character(bp$type), levels=c("non", "hub"))
#bp$hub <- factor(as.character(bp$hub), levels=unique(bp$hub))
rownames(bp) <- rownames(hmap) <- NULL

save(hmap, file=paste0(out.dir, "/", gcb, "_", out.id, "_hmap.RData"))
save(bp, file=paste0(out.dir, "/", gcb, "_", out.id, "_bp.RData"))

#-------------------------------------------------------------------------------

stat.v <- unique(hmap$stat)

for(stat in stat.v){ 
  
  #-------------------Heatmap
  
  mx <- hmap[hmap$stat==stat,]
  rwnme.hmap <- mx$h.id
  mx <- data.matrix(mx[,!colnames(mx)%in%c("h.id", "stat")])
  dimnames(mx)[[1]] <- rwnme.hmap
  mx[!is.finite(mx)] <- NA
  # Remove row/hub with all NAs
  mx <- mx[rowSums(is.na(mx))<ncol(mx),]
  # Order mx based on rowSum
  mx <- mx[order(rowSums(mx, na.rm=T), decreasing=T),]
  
  pdf(paste0(out.dir, "/", out.id, "_", stat, "_hmap.pdf"),
      width=30, height=30)
  heatmap.2(x=mx, Rowv=F, Colv=T, dendrogram="col", scale="none",
            main=paste0(out.id, "_", stat), margins=c(10,10), 
            trace="none", na.color="#4d494c", xlab=NULL, ylab=NULL, 
            breaks=c(-1,-0.5,0,0.5,1), col=c("#63a8e6", "gray90", "gray90", "#f0899a"), 
            # Note that rows can be lost when text size is set to be big
            # relative to canvas size (WEIRD)
            key=T, key.title=NA, keysize=0.5, cexRow=0.5, 
            distfun=function(x) dist(x, method="binary"))
  dev.off()
  rm(mx); gc()
  
  #-------------------Individual hub vs. boot plot per hub
  
  plot.TF <- grepl(x=names(p.lst), pattern=paste0(".", stat, "_"), fixed=T)
  p.arr <- ggarrange(plotlist=p.lst[plot.TF], nrow=5, ncol=5)
  ggexport(p.arr, height=50, width=50,
           filename=paste0(out.dir, "/", out.id, "_", stat, ".pdf"))
  rm(plot.TF, p.arr)
  
}
rm(hmap); gc()

#-------------------Boxplot with p-value from Wilcoxon paired signed test

stat.v.len <- length(stat.v)
pdf(file=paste0(out.dir, "/", out.id, "_outlier", showOutlier, "_bp.pdf"), width=40, 
    height=10*stat.v.len)
par(mfrow=c(stat.v.len,1))

HUB.len <- length(unique(bp$hub))
axis.x <- seq(from=1.5, to=HUB.len*2, by=2)
hub.TF <- bp$type=="hub"
WXTEST <- list()
stat.v <- unique(bp$stat)

for(stat in stat.v){
  
  stat.TF <- bp$stat==stat
  
  #-------------------Wilcox paired t-test (nonparametric) hub vs. nonhub
  
  wx.pval <- NULL
  trend <- NULL
  
  for( hub in unique(bp$hub) ){
    
    hub1.TF <- bp$hub==hub
    
    non.med <- median(bp[stat.TF & hub1.TF & !hub.TF, "value"], na.rm=T)
    hub.med <- median(bp[stat.TF & hub1.TF & hub.TF,"value"], na.rm=T)
    
    non.mean <- mean(bp[stat.TF & hub1.TF & !hub.TF, "value"], na.rm=T)
    hub.mean <- mean(bp[stat.TF & hub1.TF & hub.TF,"value"], na.rm=T)
    
    if( (hub.med<non.med) & (hub.mean<non.mean) ){
      trend[hub] <- -1
    } else if( (hub.med>non.med) & (hub.mean>non.mean) ){
      trend[hub] <- 1
    } else {
      # Tie or mean and median trend inconsistent
      trend[hub] <- 0
    }
      
    nh.v <- bp[stat.TF & hub1.TF & !hub.TF,"value"]
    h.v <- bp[stat.TF & hub1.TF & hub.TF,"value"]
    
    # Count ties because they should be removed to avoid this warning;
    # Also to make sure that this is not the case for most of the data
    ind.ties <- which(nh.v==h.v)
    ind.ties.len <- length(ind.ties)
    if( length(ind.ties)>0 ){
      nh.v <- nh.v[-ind.ties]
      h.v <- h.v[-ind.ties]
      print(paste0(hub, "_", stat, ": ", ind.ties.len, " tie datapoints removed."), quote=F)
    }
    
    # Just use two.sided then display non-hub and hub distribution beside each other
    # to just trend visually. Error can be there if nh.v and/or h.v are mostly 0s.
    wx.pval[hub] <- wilcox.test(x=nh.v, y=h.v, paired=T, 
                                alternative="two.sided")$p.value
   
    rm(hub1.TF, non.med, hub.med, nh.v, h.v, ind.ties.len, ind.ties)
    
    print(hub, quote=F)
    
  }
  
  WXTEST[[stat]] <- data.frame(hub=unique(bp$hub), wx.pval=wx.pval, trend=trend,
                               stringsAsFactors=F)
  sig.v <- wx.pval
  sig.v[wx.pval<0.001] <- "***"
  sig.v[wx.pval<0.01 & wx.pval>=0.001] <- "**"
  sig.v[wx.pval<0.05 & wx.pval>=0.01] <- "*"
  sig.v[wx.pval>=0.05] <- "n.s."
  
  #-------------------Plot
  
  bp$hub <- factor(x=bp$hub, levels=unique(bp$hub))
  bp$type <- factor(x=bp$type, levels=c("non", "hub"))
  
  col.x <- unname(col.v[ as.character(WXTEST[[stat]][levels(bp$hub),"trend"]) ])
  col.x <- as.vector(rbind(rep(x="gray60", times=length(col.x)), col.x))
  
  boxplot(value~type*hub, outline=showOutlier, data=bp[stat.TF,], boxwex=0.6, 
          xlab=NULL, ylab=stat, main=paste0(out.id, "_", stat, "\n",
                                            paste(sig.v, collapse="|")), 
          cex.axis=1, col=col.x, xaxt="n")
  axis(side=1, at=axis.x, labels=levels(bp$hub), cex.axis=0.5, mgp=c(3, 1.4, 0))
  # Line separating two boxplots per Cp
  abline(v=seq(from=0.5,to=HUB.len*2+1, by=2), lty=1, lwd=2, col="gray50")
  legend("topright", legend=c("non", names(col.v)), col=c("gray60", col.v), pch=15, 
         bty="n", bg="white", pt.cex=2, cex=2, horiz=F, inset=c(0,0), xpd=T)
  
  rm(sig.v, stat.TF, wx.pval, trend)
  
} # stat.v for loop end

dev.off()

rownames(WXTEST) <- NULL
save(WXTEST, file=paste0(out.dir, "/", out.id, "_wxtest.RData"))

# rm(list=ls()); gc()



