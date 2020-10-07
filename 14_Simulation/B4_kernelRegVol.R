################################################################################
# Calculate confusion matrix-derived metrics
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
csv.dir = paste0(wk.dir, "/out_compare/test")
out.dir = paste0(wk.dir, "/out_kernelRegVol")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1" 
src.id = "whole_maskMidSquare_gap50up"
out.id = "SIMvsCsnormCp_test"
confMxMetric.v = c("MCC", "FDR", "TPR", "TNR", "PPV", "NPV")

# Kernel regression estimate parameters
x.grid=seq(0,1,0.01)
y.grid=seq(0,1,0.01)
# Fixed bandwith for x and y
bws.v=c(0.05,0.05)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(np)
library(plotly)
library(RColorBrewer)
source(paste0(wk.dir, "/lib/confusionMxMetric.R"))
source(paste0(wk.dir, "/lib/doMvKernelReg.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
setwd(out.dir)

csv.v <- list.files(path=csv.dir, full.names=FALSE)
csv.len <- length(csv.v)

dta <- list()
for(c in 1:csv.len){
  
  csv <- csv.v[c]
  x <- read.csv(file=paste0(csv.dir, "/", csv), header=TRUE, stringsAsFactors=FALSE)

  # Find cut-off ranges for subject and referencec
  x$c.offsubj <- x$c.offsubj/max(x$c.offsubj)
  x$c.offref  <- x$c.offref/max(x$c.offref)
  
  # Filter based on cut-off ranges
  x <- x[x$c.offsubj>=0 & x$c.offref>=0,]
  
  dta[[c]] <- mapply(metric=confMxMetric.v, FUN=confMxMetric, MoreArgs=list(CONFMX=x))
  dta[[c]] <- cbind(subj=x$c.offsubj, ref=x$c.offref, dta[[c]])
  rm(x)
  
}
csv.v <- unlist(strsplit(x=csv.v, split=".csv", fixed=TRUE))
names(dta) <- csv.v

surf.v <- gsub(x=csv.v, pattern=paste0(gcb, "_", chr, "_"), replacement="", fixed=TRUE)
out.name <- paste(gcb, chr, out.id, sep="_")
KERREG <- list()

x.len <- length(x.grid)
y.len <- length(y.grid)
coul <- colorRampPalette(brewer.pal(11, "Spectral"))(csv.len)
for(metric in confMxMetric.v){
  
  # Kernel regression estimate
  x <- lapply(X=names(dta), FUN=function(nme){
    mx <- dta[[nme]]
    doMvKernelReg(X=mx[,1:2], Y=mx[,metric], x.grid=x.grid, 
                  y.grid=y.grid, bws.v=bws.v, 
                  plotTitle=paste0(nme, "_", metric, "_", out.id), 
                  plotPath=paste0(out.dir, "/", nme, "_", metric, "_",  
                                  out.id, "_kernelReg.html"))
  })
 
  # Plot surfaces together
  p <- plot_ly(showscale=TRUE)
  p <- layout(p, scene=list(xaxis=list(range=c(0,1)), yaxis=list(range=c(0,1))),
              title=paste0(out.name, "_", metric))
  for( i in 1:length(x) ){
    p <- add_trace(p=p, z=matrix(x[[i]][["mean"]], nrow=x.len, ncol=y.len),
                   x=x.grid, y=y.grid, type="surface", opacity=0.8,
                   colorscale=list(c(0, 1), c("white", coul[i])), 
                   colorbar=list(title=surf.v[i], len=0.2)
                   )
  }
  htmlwidgets::saveWidget(widget=as_widget(p), 
                          file=paste0(out.dir, "/", gcb, "_", chr, "_", metric, "_",  
                                      out.id, "_kernelReg.html"))
  
  # Calculate sum of estimate and save data
  y <- lapply(X=x, FUN=function(obj){
    len <- length(obj[["mean"]])
    est <- sum(obj[["mean"]])/len
    err <- sum(obj[["merr"]])/len
    c(lower.b=est-err, reg.est=est, upper.b=est+err)
  })
  y <- do.call("rbind", y)
  
  KERREG[[metric]] <- cbind.data.frame(id=csv.v, metric=metric, y, 
                                       rank.est=rank(-y[,"reg.est"], 
                                                     ties.method="average"),
                                       stringsAsFactors=FALSE)
  
}

KERREG <- do.call("rbind.data.frame", c(KERREG, stringsAsFactors=FALSE))
write.csv(KERREG, file=paste0(out.dir, "/", out.name, ".csv"),
          row.names=FALSE) 

# rm(list=ls()); gc()

