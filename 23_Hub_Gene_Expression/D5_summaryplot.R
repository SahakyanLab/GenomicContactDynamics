################################################################################
# Plot relating mean,sd, var and sd/mean to level of expression of hub genes.
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
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
bp.dir = wxtest.dir = paste0(wk.dir, "/out_hubVsNonhub/All_LTr/exprCutoff_05_acchubs")
out.dir = paste0(wk.dir, "/out_summaryplot/All_LTr/exprCutoff_05_acchubs")
### OTHER SETTINGS #############################################################
# Should match exprDataPath
src.id = "data2"

gcb = "min2Mb"
hub.id = "All_topCP3_gapBin50"

# Test statistic the grouping of hubs will be based on
stat.group = "SDdivMEAN"
# Generate plot for which test statistics:
stat.plot.v = c("MEAN", "MEDIAN", "SDdivMEAN", "NE", "LE", "ME", "HE")
col.v <- c(MEAN="#32a8a4", MEDIAN="#32a8a4", SDdivMEAN="#32a8a4", NE="gray90", 
           LE="#2171B5", ME="#9ECAE1", HE="#A50F15")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(gplots)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(gcb, "_", hub.id, "_", src.id)

# Boxplot data
load(file=paste0(bp.dir, "/", gcb, "_", out.id, "_bp.RData"))
bp <- bp[bp$type=="hub",]

# WXTEST
load(file=paste0(wxtest.dir, "/", out.id, "_wxtest.RData"))
# Indicate if the trend for each test statistic is not significant or sig.less/
# greater based on wilcoxon test.
WXTEST <- lapply(X=WXTEST[unique(c(stat.group, stat.plot.v))], FUN=function(wxstat.df){
  wxstat.df$group <- "n"
  wxstat.df$group[wxstat.df$wx.pval<0.05 & wxstat.df$trend==-1]  <- "l"
  wxstat.df$group[wxstat.df$wx.pval<0.05 & wxstat.df$trend==1] <- "g"
  wxstat.df$group <- factor(wxstat.df$group, levels=c("n", "l", "g"))
  return(wxstat.df)
})

group.df <- WXTEST$SDdivMEAN
bp <- merge(x=bp, y=group.df[,c("hub", "group")], all.x=TRUE, by="hub")

hub.v <- unique(bp$hub)

expr.level.v <- c("NE", "LE", "ME", "HE")

p.lst <- list()
for(stat in stat.plot.v){
  
  text.df <- WXTEST[[stat]]
  text.df$label <- text.df$group
  text.df$group <- group.df$group
  
  if(!stat%in%expr.level.v){
    p <- ggplot(data=bp[bp$stat==stat,], aes(x=hub, y=value)) +
      geom_boxplot(fill=col.v[stat]) +
      geom_text(aes(x=hub, label=label), y=1.5, data=text.df, size=3,
                vjust=1)+
      labs(x=NULL, y=stat, 
           title=paste0(out.id, "_", stat, "_group=", stat.group)) + 
      facet_grid(.~group, 
                 labeller=labeller(group=c(n="n.s.", l="less", g="greater"))) +
      bgr2 +
      theme(axis.text.x = element_text(size=5, angle=90, colour="black"),
            panel.grid.major.x=element_line(colour="gray90"))
    ggsave(filename=paste0(out.dir, "/", out.id, "_", stat, ".pdf"), plot=p,
           width=20, height=5); rm(p)
    p.lst <- list()
  } else {
    p.lst[[stat]] <- ggplot(data=bp[bp$stat==stat,], aes(x=hub, y=value)) +
      geom_boxplot(fill=col.v[stat]) +
      scale_y_continuous(limits=c(0,1.1), breaks=seq(0,1,0.2)) +
      geom_text(aes(x=hub, label=label), y=1.1, data=text.df, size=3,
                vjust=1)+
      labs(x=NULL, y="Gene fraction", 
           title=paste0(out.id, "_", stat, "_group=", stat.group)) + 
      facet_grid(.~group, 
                 labeller=labeller(group=c(n="n.s.", l="less", g="greater"))) +
      bgr2 +
      theme(axis.text.x = element_text(size=5, angle=90, colour="black"),
            panel.grid.major.x=element_line(colour="gray90"))
  }
  
  print(stat)
  
}

p.arr <- ggarrange(plotlist=p.lst, nrow=4, ncol=1)
ggexport(p.arr, height=20, width=20,
         filename=paste0(out.dir, "/", out.id, ".pdf"))

# rm(list=ls()); gc()

  
