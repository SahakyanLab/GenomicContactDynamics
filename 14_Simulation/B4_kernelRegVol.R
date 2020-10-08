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
confMxMetric.v = "MCC" #c("MCC", "FDR", "TPR", "TNR", "PPV", "NPV")

# Cut-off range to consider
s.range = 'subj.range = c(-0.0001,0.004)'
r.range = 'ref.range = c(0,5)' 

# Kernel regression estimate parameters
step = 'step.v=c(x.grid=0.01, y.grid=0.01)'
bws = 'bws.v=c(0.05, 0.05)' # Fixed bandwith for x and y

# Bootstrapping of grid estimates based on kernel regression estimate error
NBOOT=10
SEED=3412
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(np)
library(plotly)
library(RColorBrewer)
library(reshape2) # For acast() in getVUS.R
source(paste0(wk.dir, "/lib/confusionMxMetric.R"))
source(paste0(wk.dir, "/lib/doMvKernelReg.R"))
source(paste0(wk.dir, "/lib/drawGridEst.R"))
source(paste0(wk.dir, "/lib/getVUS.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
setwd(out.dir)

#print(paste(gcb, chr, src.id, out.id, sep="_"), quote=FALSE)
out.name <- paste0(paste(gcb, chr, src.id, out.id, sep="_"), "_seed", SEED, "_nboot", NBOOT)
print(paste0(out.name, "..."), quote=FALSE)
print(step, quote=FALSE)
print(bws, quote=FALSE)
print(s.range, quote=FALSE)
print(r.range, quote=FALSE)

eval(parse(text=step))
eval(parse(text=bws))
eval(parse(text=s.range))
eval(parse(text=r.range))

csv.v <- list.files(path=csv.dir, full.names=FALSE)
csv.len <- length(csv.v)
dta <- list()
for(c in 1:csv.len){
  
  csv <- csv.v[c]
  x <- read.csv(file=paste0(csv.dir, "/", csv), header=TRUE, stringsAsFactors=FALSE)

  # Filter based on subj and ref cut-off ranges
  if( is.null(subj.range ) ){
    subj.range <- c( min(x$c.offsubj), max(x$c.offsubj) )
  }
  if( is.null(ref.range) ){
    ref.range <- c( min(x$c.offref), max(x$c.offref) )
  }
  out.range <- (x$c.offsubj<subj.range[1] | x$c.offsubj>subj.range[2]) |
    (x$c.offref<ref.range[1] | x$c.offref>ref.range[2])
  x <- x[out.range==FALSE,]
  rm(out.range)
  
  # Min-max normalisation of cut-off ranges -> [0,1]
  x$c.offsubj <- (x$c.offsubj-min(x$c.offsubj))/(max(x$c.offsubj)-min(x$c.offsubj))
  x$c.offref  <- (x$c.offref-min(x$c.offref))/(max(x$c.offref)-min(x$c.offref))
  
  # Calculate scores
  dta[[c]] <- mapply(metric=confMxMetric.v, FUN=confMxMetric, MoreArgs=list(CONFMX=x))
  dta[[c]] <- cbind(subj=x$c.offsubj, ref=x$c.offref, dta[[c]])
  rm(x)
  
}
csv.v <- unlist(strsplit(x=csv.v, split=".csv", fixed=TRUE))
names(dta) <- csv.v

#out.name <- paste0(paste(gcb, chr, out.id, sep="_"), "_seed", SEED, "_nboot", NBOOT)

# Plot parameters
coul <- colorRampPalette(brewer.pal(11, "Spectral"))(csv.len)
surf.v <- gsub(x=csv.v, pattern=paste0(gcb, "_", chr, "_"), replacement="", fixed=TRUE)

# Make common grid
x.grid <- seq(0, 1, step.v["x.grid"])
y.grid <- seq(0, 1, step.v["y.grid"])
x.len <- length(x.grid)
y.len <- length(y.grid)
grid.mx <- as.matrix(expand.grid(x.grid, y.grid))

KERREG <- list()
for(metric in confMxMetric.v){
  
  # Get kernel regression estimate
  KER <- lapply(X=names(dta), FUN=function(nme){
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
    drawGridEst(est.v=obj$mean, err.v=obj$merr, nboot=NBOOT, SEED=SEED)
  })
  rm(KER)
  
  # Approximate VUS for all NBOOT samples
  VUS <- lapply(X=BOOT, FUN=function(boot){
    x <- lapply(X=boot, FUN=function(est){
      getVUS( X=grid.mx, Y=est, step.v=unname(step.v[c("x.grid", "y.grid")]) )
    })
    unlist(x)
  })
  rm(BOOT)
 
  meanVUS <- unlist(lapply(X=VUS, FUN=mean))
  sdVUS <- unlist(lapply(X=VUS, FUN=sd))
  KERREG[[metric]] <- cbind.data.frame(id=csv.v, metric=metric, nboot=NBOOT, 
                                       mean=meanVUS, sd=sdVUS,
                                       rank.est=rank(-meanVUS, 
                                                     ties.method="average"),
                                       stringsAsFactors=FALSE)
  rm(meanVUS, sdVUS, VUS)
  
}

KERREG <- do.call("rbind.data.frame", c(KERREG, stringsAsFactors=FALSE))
write.csv(KERREG, file=paste0(out.dir, "/", out.name, ".csv"), row.names=FALSE) 

# rm(list=ls()); gc()

