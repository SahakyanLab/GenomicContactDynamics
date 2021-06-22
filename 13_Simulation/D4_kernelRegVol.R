################################################################################
# Convert tile plot from B2 into a 3D surface via kernel regression to compare
# different subj-ref pairs. Estimate volume of surface to serve as performance 
# measure (similar to AUC of ROC) via trapezoidal rule for double integration. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
csv.dir = paste0(wk.dir, "/out_compare_sd_csv")
out.dir = paste0(wk.dir, "/out_kernelRegVol")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1" 
src.id = "whole_maskMidSquare_gap50up"
out.id = "SIMvsCsnorm"
confMxMetric.v = "MCC" #c("MCC", "FDR", "TPR", "TNR", "PPV", "NPV")

# Cut-off range to consider
s.range = 'subj.range = c(-0.0001,0.004)'
r.range = 'ref.range = c(0,5)' 

# Kernel regression estimate parameters
step = 'step.v=c(xgrid=0.01, ygrid=0.01)'
bws = 'bws.v=c(subj=0.05, ref=0.05)' # Fixed bandwidth for each dimension

# Bootstrapping of grid estimates based on kernel regression estimate error
NBOOT = 10
nCPU = 1
SEED = 902

# Get optimal bandwidth only?
getBWonly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(np) # For kernel regression, npregbw() and npreg()
library(plotly)
library(RColorBrewer)
library(reshape2) # For acast() in getVUS.R
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(wk.dir, "/lib/confusionMxMetric.R"))
source(paste0(wk.dir, "/lib/doMvKernelReg.R"))
source(paste0(wk.dir, "/lib/drawGridEst.R"))
source(paste0(wk.dir, "/lib/getVUS.R"))
### FUNCTION ###################################################################
getminBW <- function(v, v.grid){
  v <- sort(v)
  v.grid <- sort(v.grid)
  max(
    sapply(X=v, simplify=TRUE, FUN=function(x){
      min( abs(x-v.grid) )
    })
  )
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# For plotly
setwd(out.dir)

print(step, quote=FALSE)
print(bws, quote=FALSE)
print(s.range, quote=FALSE)
print(r.range, quote=FALSE)

eval(parse(text=step))
eval(parse(text=bws))
eval(parse(text=s.range))
eval(parse(text=r.range))

step <- paste( paste(names(step.v), step.v, sep=""), collapse="_") 
param.id <- paste0("seed", SEED, "_nboot", NBOOT, "_step", step)
out.name <- paste0(paste(gcb, chr, src.id, out.id, sep="_"), "_", param.id); rm(step)
print(paste0(out.name, "..."), quote=FALSE)

# Make common grid for kernel regression estimation of score and VUS calculation
x.grid <- seq(0, 1, step.v["xgrid"])
y.grid <- seq(0, 1, step.v["ygrid"])
step.v <- unname(step.v[c("xgrid", "ygrid")])
x.len <- length(x.grid)
y.len <- length(y.grid)
grid.mx <- as.matrix(expand.grid(x.grid, y.grid))

csv.v <- list.files(path=csv.dir, full.names=FALSE)
csv.v <- unlist(strsplit(x=csv.v, split=".csv", fixed=TRUE))
csv.len <- length(csv.v)
dta <- BW <- list()
for(c in 1:csv.len){
  
  x <- read.csv(file=paste0(csv.dir, "/", csv.v[c], ".csv"), header=TRUE, 
                stringsAsFactors=FALSE)

  # Filter based on subj and ref cut-off ranges
  if( is.null(subj.range ) ){
    subj.range <- c( min(x$c.offsubj), max(x$c.offsubj) )
  }
  if( is.null(ref.range) ){
    ref.range <- c( min(x$c.offref), max(x$c.offref) )
  }
  out.range <- (x$c.offsubj<subj.range[1] | x$c.offsubj>subj.range[2]) |
    (x$c.offref<ref.range[1] | x$c.offref>ref.range[2])
  x <- x[!out.range,]
  rm(out.range)
  
  # Min-max normalisation of cut-off ranges -> [0,1]
  x$c.offsubj <- (x$c.offsubj-min(x$c.offsubj))/(max(x$c.offsubj)-min(x$c.offsubj))
  x$c.offref  <- (x$c.offref-min(x$c.offref))/(max(x$c.offref)-min(x$c.offref))
  
  # Calculate scores
  dta[[c]] <- mapply(metric=confMxMetric.v, FUN=confMxMetric, MoreArgs=list(CONFMX=x))
  dta[[c]] <- cbind(subj=x$c.offsubj, ref=x$c.offref, dta[[c]])
  rm(x)
  
  # Find and store optimal bandwidth per metric, per matrix pair
  # Optimal such that kernel function covers at least one neighbor
  BW[[c]] <- sapply(X=confMxMetric.v, simplify=TRUE, FUN=function(metric){
    x <- dta[[c]][is.finite(dta[[c]][,metric]),1:2]
    sBW <- getminBW(v=unique(x[,"subj"]), v.grid=x.grid)
    rBW <- getminBW(v=unique(x[,"ref"]), v.grid=y.grid)
    return( c(minBWsubj=sBW, minBWref=rBW) )
  })
  BW[[c]] <- cbind.data.frame(id=csv.v[c], metric=colnames(BW[[c]]), t(BW[[c]]))

}
gc()
names(dta) <- csv.v

write.csv(do.call("rbind.data.frame", BW), 
          file=paste0(out.dir, "/", out.name, "_minBW.csv"), row.names=FALSE)
rm(BW)

if(getBWonly){
  
  print("Only determine optimal bandwidth.", quote=FALSE)
  
} else {
  
  bws <- paste( paste(names(bws.v), bws.v, sep=""), collapse="_") 
  param.id <- paste0(param.id, "_bws", bws)
  out.name <- paste0(out.name, "_bws", bws); rm(bws)
  
  # Plot parameters
  coul <- colorRampPalette(brewer.pal(11, "Spectral"))(csv.len)
  surf.v <- gsub(x=csv.v, pattern=paste0(gcb, "_", chr, "_"), replacement="", fixed=TRUE)
  
  KERREG <- list()
  for(metric in confMxMetric.v){
    
    # Get kernel regression estimate
    KER <- lapply(X=names(dta), FUN=function(nme){
      mx <- dta[[nme]]
      doMvKernelReg(X=mx[,1:2], Y=mx[,metric], x.grid=x.grid, 
                    y.grid=y.grid, bws.v=unname(bws.v[c("subj", "ref")]), 
                    plotTitle=paste0(nme, "_", param.id, "_", metric, "_", out.id), 
                    plotPath=paste0(out.dir, "/", nme, "_", param.id, "_", 
                                    metric, "_", out.id, "_kernRegSurf.html"))
    })
    
    # Plot surfaces together
    p <- plot_ly(showscale=TRUE)
    p <- layout(p, scene=list(xaxis=list(range=c(0,1)), yaxis=list(range=c(0,1))),
                title=paste0(out.name, "_", metric))
    for( i in 1:length(KER) ){
      p <- add_trace(p=p, z=matrix(KER[[i]][["mean"]], nrow=x.len, ncol=y.len),
                     x=x.grid, y=y.grid, type="surface", opacity=0.8,
                     colorscale=list(c(0, 1), c("white", coul[i])), 
                     colorbar=list(title=surf.v[i], len=0.2)
      )
    }
    htmlwidgets::saveWidget(widget=as_widget(p), 
                            file=paste0(out.dir, "/", out.name, "_", metric, 
                                        "_surfCombined.html"))
    rm(p)
    
    # Get NBOOT samples of kernel estimates at common grid
    BOOT <- lapply(X=KER, FUN=function(obj){
      drawGridEst(est.v=obj$mean, err.v=obj$merr, nboot=NBOOT, SEED=SEED, nCPU=nCPU)
    })
    rm(KER)
    
    # Approximate VUS for all NBOOT samples
    VUS <- lapply(X=BOOT, FUN=function(boot.mx){
      x <- apply(X=boot.mx, MARGIN=2, FUN=function(est){
        getVUS(X=grid.mx, Y=est, step.v=step.v)
      })
      return(x)
    })
    rm(BOOT)
    VUS <- do.call("cbind", VUS)
    
    meanVUS <- colMeans(x=VUS)
    sdVUS <- apply(X=VUS, MARGIN=2, FUN=sd)
    KERREG[[metric]] <- cbind.data.frame(id=csv.v, metric=metric, nboot=NBOOT, 
                                         mean=meanVUS, sd=sdVUS,
                                         rank.est=rank(-meanVUS, 
                                                       ties.method="average"),
                                         stringsAsFactors=FALSE)
    rm(meanVUS, sdVUS, VUS); gc()
    
  }
  
  KERREG <- do.call("rbind.data.frame", c(KERREG, stringsAsFactors=FALSE))
  write.csv(KERREG, file=paste0(out.dir, "/", out.name, "_VUS.csv"), row.names=FALSE) 
  
}

# rm(list=ls()); gc()

