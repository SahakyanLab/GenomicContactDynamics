################################################################################
# Multivariate kernel regression
# Source: https://bookdown.org/egarpor/NP-UC3M/kre-ii-multmix.html
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(np)
# library(plotly)
### FUNCTION ###################################################################
doMvKernelReg <- function(X, Y, x.grid, y.grid, plotTitle, bws.v=c(0.01,0.25),
                          plotPath=paste0(out.dir, "/", out.name, ".html")
                          ){
  
  x.len <- length(x.grid)
  y.len <- length(y.grid)
  grid.mx <- as.matrix(expand.grid(x.grid, y.grid))
  dimnames(grid.mx)[[2]] <- c("x", "y")
  
  BW <- np::npregbw(xdat=X, ydat=Y, bws=bws.v, bwtype="fixed", regtype="lc",
                    # bwscaling=FALSE means bws are interpreted as raw bandwidths
                    bwscaling=FALSE, bandwidth.compute=FALSE, ckertype="gaussian",
                    na.action=na.omit)
  print(paste0("Bandwidths: ", BW$bw), quote=FALSE)
  NPREG <- np::npreg(bws=BW, exdat=grid.mx, gradients=TRUE)
  
  # Save individual surface plots
  p <- plot_ly(x=x.grid, y=y.grid, z=matrix(NPREG$mean, nrow=x.len, ncol=y.len))
  p <- add_surface(p=p)
  p <- layout(p, scene=list(xaxis=list(range=c(0,1)), yaxis=list(range=c(0,1))),
              title=plotTitle)
  htmlwidgets::saveWidget(widget=as_widget(p), file=plotPath)
  
  return(NPREG)
  
}
################################################################################
