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
  
  png(filename=paste0(out.dir, "/", out.name, "_VsMaxCp.png"), 
      height=300*10, width=300*10, res=300)
  
  boxplot(formula=corval~maxCp, data=CORCP.DF, outline=F, xlab="Cp", ylab=cor.meth, 
          main=paste0(out.name, TEST$test.id), col="#FDC776", cex.main=0.8, ylim=c(-1,1))
  
  dev.off()
  
}

# rm(list=ls()); gc()