################################################################################
# Boxplot for geneExprVsCp.R
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(reshape2)
# library(ggplot2)
# source(paste0(lib, "/GG_bgr.R"))
# source(paste0(lib, "/compareTwoDist.R"))
### FUNCTION ###################################################################
plotgeneExprVsCp <- function(CPGNEXPR.DF = 'CPGNEXPR.DF',
                             manner='acrossTissues | perTissue',
                             out.dir = 'directory',
                             out.name = 'filename',
                             plottype = 'boxplot | foldchangeplot',
                             plottitle = 'test',
                             statToPlot = 'c("nNE", "nHE")',
                             
                             funxforfoldchange = 'mean | median',
                             box.outliers = 'T/F',
                             dim.boxplot = 'for par(mfrow)'){
  
  # Specify colours for stats
  #col.v = c(MEAN="#d9c741", MEDIAN="#edb340", SD="#32a864", SDdivMEAN="#32a8a4", 
  #          nND="gray50", nNE="gray70", nLE="#2171B5", nME="#9ECAE1", nHE="#A50F15")
  col.v = c(MEAN="#2171B5", MEDIAN="#edb340", SD="#32a864", SDdivMEAN="#b592d4", 
            nND="gray50", nNE="gray70", nLE="#2171B5", nME="#9ECAE1", nHE="#A50F15")
  #"#a175c9"
  if( any(!statToPlot%in%colnames(CPGNEXPR.DF)) ){
    stop("plotgeneExprVsCp(): Invalid statToPlot element.")
  } else {
    CPGNEXPR.DF <- CPGNEXPR.DF[,unique(c("Cp", "id", statToPlot))]
  }
  
  CPGNEXPR.DF <- reshape2::melt(data=CPGNEXPR.DF, id=c("Cp", "id"))
  
  # Fix types of columns
  Cp.v <- unique(as.character(setdiff(CPGNEXPR.DF$Cp, "HiC_all")))
  Cp.v <- c( "HiC_all", sort(as.numeric(Cp.v), decreasing=F))
  CPGNEXPR.DF$Cp <- factor(x=as.character(CPGNEXPR.DF$Cp), levels=Cp.v)
  CPGNEXPR.DF$id <- as.character(CPGNEXPR.DF$id)
  CPGNEXPR.DF$variable <- as.character(CPGNEXPR.DF$variable)

  if(plottype=="boxplot"){
    
    pdf(file=paste0(out.dir, "/", out.name, "_outlier", box.outliers, "_", 
                    plottype, ".pdf"), height=dim.boxplot[[1]]*10, width=dim.boxplot[[2]]*10)
    par(mfrow=dim.boxplot)
    
    for(stat in statToPlot){
      
      if( stat%in%c("nND", "nNE", "nLE", "nME", "nHE") ){
        ylim.val <- c(0,1)
      } else {
        ylim.val <- NULL
      }
      
      stat.TF <- CPGNEXPR.DF$variable==stat
      TEST <- compareTwoDist(x=CPGNEXPR.DF$value[stat.TF & as.character(CPGNEXPR.DF$Cp)%in%c("1", "2", "3")],
                             y=CPGNEXPR.DF$value[stat.TF & as.character(CPGNEXPR.DF$Cp)%in%c("19", "20", "21")])
      
      boxplot(value~Cp, data=CPGNEXPR.DF[stat.TF,],
              outline=box.outliers, boxwex=0.6, xlab=expression("c"["p"]), 
              ylab=stat, cex.axis=1.5, main=paste0(TEST$test.id, "_x=Cp1To3,y=Cp19To21"), 
              xaxt="n", col=col.v[[stat]], ylim=ylim.val)
      mtext(text=plottitle, side=3, line=2, cex=1)
      axis(side=1, at=1:length(levels(CPGNEXPR.DF$Cp)), labels=levels(CPGNEXPR.DF$Cp), 
           cex=1.5)
      rm(stat.TF)
      
    }
    
    dev.off()
    
  } else if(plottype=="foldchangeplot"){
    
    eval(parse(text=paste0(
      
      'df <- aggregate(x=CPGNEXPR.DF$value, by=list(CPGNEXPR.DF$variable, 
             as.character(CPGNEXPR.DF$Cp)), FUN=', funxforfoldchange, ', na.rm=T)'
      
    )))
    rm(CPGNEXPR.DF)
    colnames(df) <- c("stat", "Cp", "value")
    
    df$FC <- log2(df$value/df$value[as.character(df$Cp)=="HiC_all"])
    
    df <- df[df$Cp!="HiC_all",] 
    Cp.v <- setdiff(Cp.v, "HiC_all")
    
    df$Cp <- factor(x=df$Cp, levels=Cp.v)
    df$stat <- as.character(df$stat)
    df$stat <- factor(x=df$stat, levels=union(statToPlot, unique(df$stat)))
    
    p <- ggplot(data=df, aes(x=Cp, y=FC, group=stat, colour=stat)) +
      geom_line(size=2) +
      geom_point(size=2.5) +
      scale_y_continuous(breaks=seq(from=-2, to=0.5, by=0.5), limits=c(-1,0.5)) + 
      scale_colour_manual(values=col.v[levels(df$stat)]) + 
      labs(y=paste0("fc ", funxforfoldchange, " wrt HiC_all"), colour=NULL, 
           title=plottitle) +
      bgr2 +
      theme(axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10))
    
    ggsave(filename=paste0(out.dir, "/", out.name, "_funx", funxforfoldchange,
                           "_", plottype, ".pdf"), 
           height=10, width=10, plot=p)
    
  } else {
    stop("plotgeneExprVsCp(): Invalid plottype argumet.")
  }
  
}

# rm(list=ls()); gc()