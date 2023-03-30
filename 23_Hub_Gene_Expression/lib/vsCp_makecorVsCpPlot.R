################################################################################
# Function to generate correlation coefficient vs. Cp dimension
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(lib, "/compareTwoDist.R"))
### FUNCTION ###################################################################
makecorVsCpPlot <- function(CORCP.DF, CpsToCompare, out.dir, out.name){
  
  CORCP.DF <- CORCP.DF[!is.na(CORCP.DF$maxCp),]
  
  TEST <- compareTwoDist(x=CORCP.DF$corval[ CORCP.DF$maxCp %in% CpsToCompare[[1]] ],
                         y=CORCP.DF$corval[ CORCP.DF$maxCp %in% CpsToCompare[[2]] ])
  
  Cp.v <- as.character(sort(unique(CORCP.DF$maxCp)))
  CORCP.DF$maxCp <- factor(x=as.character(CORCP.DF$maxCp), levels=Cp.v)
  
  med.val <- median(CORCP.DF$corval[as.character(CORCP.DF$maxCp) == "1"])
  
  p <- ggplot(data=CORCP.DF, aes(x=maxCp, y=corval)) +
    geom_hline(yintercept=med.val, linetype="solid", colour="gray50") + 
    stat_boxplot(geom="errorbar", width=0.5, outlier.colour=adjustcolor("black", 0.1)) + 
    geom_boxplot(fill=adjustcolor("#FDC776", alpha=1)) +
    labs(title=paste0(out.name, TEST$test.id)) + 
    scale_y_continuous(limits=c(-1,1)) +
    bgr2 
  
  ggsave(filename=paste0(out.dir, "/", out.name, "_VsMaxCp.png"),
         width=10, height=10)
  
  #png(filename=paste0(out.dir, "/", out.name, "_VsMaxCp.png"), 
  #    height=300*10, width=300*10, res=300)
  
  #boxplot(formula=corval~maxCp, data=CORCP.DF, outline=F, xlab="Cp", ylab=cor.meth, 
  #        main=paste0(out.name, TEST$test.id), col="honeydew3", cex.main=0.8, ylim=c(-1,1))
  #abline(h=0.5)
  
  #dev.off()
  
  #p <- ggplot(data=CORCP.DF, aes(x=maxCp, y=corval)) +
  #  geom_violin(scale="width", col="darkblue", alpha=0.5, trim=T, lwd=2) +
  #  bgr2
  
}

# rm(list=ls()); gc()