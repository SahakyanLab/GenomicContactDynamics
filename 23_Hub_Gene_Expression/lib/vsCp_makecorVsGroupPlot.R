################################################################################
# Function to plot correlation coefficient distributions of dynamic/persistent 
# max Cp and non-hub/hub gene pairs
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(ggplot2)
# source(paste0(wk.dir, "/lib/compareManyDist.R"))
### FUNCTION ###################################################################
makecorVsGroupPlot <- function(CORCP.DF, dynamic.Cp, persistent.Cp, 
                               out.dir, out.name){
  
  CORCP.DF <- CORCP.DF[ !is.na(CORCP.DF$hubTF) & !is.na(CORCP.DF$maxCp),]
  #CORCP.DF$maxCp <- as.numeric(as.character(CORCP.DF$maxCp))
  
  # Apply Cp range for dynamic and persistent max Cp
  CORCP.DF <- CORCP.DF[CORCP.DF$maxCp %in% c(dynamic.Cp, persistent.Cp),]
  
  CORCP.DF$group <- NA  
  CORCP.DF$group[CORCP.DF$maxCp %in% dynamic.Cp] <- paste(range(dynamic.Cp), collapse="-")
  CORCP.DF$group[CORCP.DF$maxCp %in% persistent.Cp] <- paste(range(persistent.Cp), collapse="-")
  CORCP.DF$Cpgroup <- CORCP.DF$group
  
  CORCP.DF <- CORCP.DF[order(CORCP.DF$hubTF), ]
  CORCP.DF <- CORCP.DF[order(CORCP.DF$Cpgroup, decreasing=T), ]
  
  CORCP.DF$group <- paste0("Cp", CORCP.DF$Cpgroup, "_", CORCP.DF$hubTF, "hub")
  
  CORCP.DF$Cpgroup <- factor(as.character(CORCP.DF$Cpgroup), 
                             levels=as.character(unique(CORCP.DF$Cpgroup)))
  CORCP.DF$group <- factor(as.character(CORCP.DF$group),
                           levels=as.character(unique(CORCP.DF$group)))
    
  if( any(is.na(CORCP.DF$group)) ){
    stop("makecorVsGroupPlot(): Missing group.")
  } else {
    
    compareManyDist(xval=CORCP.DF$corval, grp=CORCP.DF$group, out.dir=out.dir, 
                    out.name=paste0(out.name, "_corval"))
    
    cols <- colorRampPalette(brewer.pal(n=11, name="Spectral"))(21)
    cols <- c(adjustcolor(cols[[1]], 1), adjustcolor(cols[[1]], 1), 
              adjustcolor(cols[[21]], 1), adjustcolor(cols[[21]], 1))
    names(cols) <- levels(CORCP.DF$group)
    
    cols.fill <- cols
    cols.fill[c(1,3)] <- "gray80"
    
    # Boxplot
    png(filename=paste0(out.dir, "/", out.name, "_corval_dynamic-persistentAndNonhub0-Hub1_box.png"), 
        height=300*10, width=300*10, res=300)
    
    boxplot(formula=corval~group, data=CORCP.DF[!is.na(CORCP.DF$maxCp),], outline=F, 
            xlab=NULL, ylab=cor.meth, main=out.name, cex.main=0.8, col=cols, 
            ylim=c(-1,1))
    
    df.mean <- stack(by(data=CORCP.DF$corval, INDICES=CORCP.DF$group, FUN=mean, na.rm=F))
    points(x=df.mean$ind, y=df.mean$values, cex=2, pch=7)
    rm(df.mean)
    
    dev.off()
    
    # Violin plot
    p <- makeViolinPlot(df=CORCP.DF, x.nme="group", y.nme="corval", sd.mult=1, col.nme="group", 
                        line.cols=cols, fill.nme="group", fill.cols=cols.fill, fill.legend="none",  
                        plot.title=out.name, ylim.val=c(-1,1), showOutlier=F, addmean=F, 
                        geom.viol.scale="area") # "area" is default, trim=T
    
    #p <- ggplot(data=CORCP.DF, aes(x=group, y=corval)) +
    #  geom_violin(scale="width", col="darkblue", alpha=0.5, trim=T, lwd=2) +
    #  bgr2
    
    ggsave(filename=paste0(out.dir, "/", out.name, "_corval_dynamic-persistentAndNonhub0-Hub1_violin.pdf"),
           plot=p, height=10, width=10)
    
    
    cols.fill <- cols
    cols.fill[c(1,3)] <- "gray80"
    # Density plot
    p <- ggplot(data=CORCP.DF, aes(x=corval)) +
      geom_density(aes(fill=group, col=group), size=1.5, alpha=0.7) + 
      scale_x_continuous(limits=c(-1,1)) + 
      scale_color_manual(values=cols) + 
      scale_fill_manual(values=cols.fill) +
      labs(title=out.name) +
      #facet_grid(Cpgroup~.) +
      bgr2
    
    ggsave(filename=paste0(out.dir, "/", out.name, "_corval_dynamic-persistentAndNonhub0-Hub1_density.pdf"),
           plot=p, height=10, width=10)
      
    
  }
  
  return(CORCP.DF)
  
}

# rm(list=ls()); gc()