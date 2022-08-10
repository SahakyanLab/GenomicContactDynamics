################################################################################
# Boxplot of calculated statistic per hub with significance from bootstrapping 
# indicated.
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
id = "exprCutoff_05_acchubs" # "exprCutoff_05_allhubs" | "exprCutoff_05_acchubs" 
# Source data and significance test result
bp.dir = wxtest.dir = paste0(wk.dir, "/out_hubVsNonhub/All_LTr/", id)
out.dir = paste0(wk.dir, "/out_summaryplot_stat/All_LTr/", id)
### OTHER SETTINGS #############################################################
src.id.v = c("data1", "data2")
gcb = "min2Mb"
hub.id = "All_topCP3_gapBin50" # Same results with All
stat = "SDdivMEAN"
# sig0trend, trend is 0 (meaning mean and median don't agree but pvalue is sig)
col.v = c(l="mediumpurple3", n="tomato", g="#51C56AFF", sig0trend="lightblue")
alpha = 0.05
outlier = F
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out0.id <- paste0(gcb, "_", hub.id)

BP <- list()
WX <- list()
for(src.id in src.id.v){
  
  out.id <- paste0(out0.id, "_", src.id)
  
  # Source
  load(file=paste0(bp.dir, "/", gcb, "_", out.id, "_bp.RData"))
  bp$type <- as.character(bp$type)
  bp$stat <- as.character(bp$stat)
  bp <- bp[bp$type=="hub" & bp$stat==stat,]
  bp$hub <- as.character(bp$hub)
  bp$hub <- factor(x=bp$hub, levels=sort(unique(bp$hub), decreasing=F))
  
  HUB.len <- length(levels(bp$hub))
  
  # WXTEST
  load(file=paste0(wxtest.dir, "/", out.id, "_wxtest.RData"))
  # Indicate if the trend for each test statistic is not significant or sig.less/
  # greater based on wilcoxon test.
  WXTEST <- WXTEST[[stat]]
  WXTEST$group <- "n"
  WXTEST$group[ WXTEST$wx.pval<alpha & WXTEST$trend==-1 ] <- "l"
  WXTEST$group[ WXTEST$wx.pval<alpha & WXTEST$trend==1  ] <- "g"
  WXTEST$group[ WXTEST$wx.pval<alpha & WXTEST$trend==0  ] <- "sig0trend"
  WXTEST$group <- factor(x=WXTEST$group, levels=c("n", "l", "g", "sig0trend"))
  rownames(WXTEST) <- WXTEST$hub
  
  bp <- merge(x=bp, y=WXTEST[,c("hub", "group", "trend")], all.x=TRUE, by="hub")
  
  # Per src.id boxplot
  pdf(file=paste0(out.dir, "/", out.id, "_", stat, "_summarybp.pdf"), 
      width=40, height=10)
  
  boxplot(value~hub, outline=outlier, data=bp, boxwex=0.6, xlab=NULL, ylab=stat, 
          main=paste0(out.id, "_alpha", alpha), cex.axis=2, xaxt="n",
          col=col.v[ as.character(WXTEST[levels(bp$hub), "group"]) ])
  axis(side=1, at=seq(1, HUB.len, 1), labels=levels(bp$hub),
       cex.axis=0.8, mgp=c(3, 1.4, 0))
  abline(v=seq(0.5, HUB.len+0.5, by=1), lty=1, col="gray50") 
  legend(x="topright", legend=names(col.v), col=col.v, lty=1, lwd=5, cex=0.8)
  
  dev.off()
  
  BP[[src.id]] <- cbind(src.id=src.id, bp)
  WX[[src.id]] <- cbind(src.id=src.id, WXTEST)
  rm(WXTEST, bp); gc()
  
} # src.id.v

BP <- do.call("rbind", BP)
rownames(BP) <- NULL

BP$hub <- factor(x=as.character(BP$hub), 
                 levels=sort(unique(as.character(BP$hub)), decreasing=F))
BP$src.id <- factor(x=as.character(BP$src.id),   
                    levels=sort(src.id.v, decreasing=F))
BP$group <- as.character(BP$group)

WX <- do.call("rbind", WX)
rownames(WX) <- paste0(as.character(WX$hub), "_", as.character(WX$src.id))
WX$group <- as.character(WX$group)

# Colours for combined boxplot
rw.id <- paste0(rep( x=levels(BP$hub), each=length(levels(BP$src.id)) ), "_",
                rep( x=levels(BP$src.id), times=length(levels(BP$hub))) 
                )
cbp.v <- col.v[ WX[rw.id, "group"] ]
trend.v <- rep(NA, length(rw.id))
trend.v[as.numeric(WX[rw.id, "trend"])==1] <- "g"
trend.v[as.numeric(WX[rw.id, "trend"])==-1] <- "l"

BP$hub <- as.character(BP$hub)
BP$hub <- gsub(x=BP$hub, pattern="_", replacement="\n", fixed=TRUE)
BP$hub <- factor(x=as.character(BP$hub), 
                 levels=sort(unique(as.character(BP$hub)), decreasing=F))

# Boxplot combining src.id
pdf(file=paste0(out.dir, "/", out0.id, "_", stat, "_summarybp.pdf"), 
    width=40, height=10)
par(mar = c(6.1, 6.1, 4.1, 4.1))

boxplot(value~src.id*hub, outline=outlier, data=BP, boxwex=0.6, xlab=NULL, 
        ylab=NULL, main=paste0(out0.id, "_alpha", alpha, "_", stat, "_order_src=",
                               paste(levels(BP$src.id), collapse=",")), 
        cex.axis=3, xaxt="n", yaxt="n", col=cbp.v, mgp=c(8, 4, 0))
axis(side=1, at=seq(1.5, HUB.len*2, 2), labels=levels(BP$hub), cex.axis=3, 
     mgp=c(8, 4.5, 0), lwd.ticks=0)
axis(side=2, las=2, cex.axis=3)
text(labels=trend.v, x=seq(1, HUB.len*2, by=1), cex=4, 
     # Using all hubs, introduces NA. Selected hubs have no NA.
     y=max(BP$value, na.rm=T))
abline(v=seq(0.5, HUB.len*2+0.5, by=2), lty=1, lwd=2, col="gray50") 
legend(x="topright", legend=names(col.v), col=col.v, lty=1, lwd=5, 
       cex=1, bty="n")

dev.off()

# rm(list=ls()); gc()

