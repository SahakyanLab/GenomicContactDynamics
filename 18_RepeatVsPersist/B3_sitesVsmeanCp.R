################################################################################
#  Plot repeat site count vs. average Cp (wmeanCp or wmeanCp0) per bin. Perform
# correlation of the two per repeat subfamily/family. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")

rep.group = "subfam"
binrep.dir = paste0(wk.dir, "/out_RepeatOverlapPerBinALL/", rep.group)
meanCp.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binWeightedMeanCp/Cp1To21")
out.dir = paste0(wk.dir, "/out_sitesVsmeanCp")

repeat.file = paste0(wk.dir, "/Repeat_rankingbyAge/repsubfamALL.csv")
agerank.file = paste0(wk.dir, "/Repeat_rankingbyAge/repsubfam.csv")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
nCPU = 1L 
cuts = 3
inclbins0sites = T
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
library(hexbin)
library(ggplot2)
library(Hmisc)
library(ggpubr)
library(viridis)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/doCorTest.R"))
source(paste0(lib, "/makeHexbinggplot.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)

repeat.v <- read.csv(file=repeat.file, header=T, stringsAsFactors=F)[,"repName"]

toExport <- c("chr.v", "meanCp.dir", "binrep.dir", "gcb", "repeat.v")
BINREPMEANCP.MX <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                           .inorder=F, .combine="rbind",
                           .export=toExport, .noexport=ls()[!ls()%in%toExport]
                        
) %op% {
  
  chunk <- sapply(X=itr, simplify=F, FUN=function(i){
    
    chr <- chr.v[i]
    load(file=paste0(meanCp.dir, "/", gcb, "_", chr, "_weightedMeanCp.RData"))
    load(file=paste0(binrep.dir,"/", chr, "_BinRep.RData"))
    
    if( identical(as.numeric(BINREP.MX[,"bins"]),
                  as.numeric(BINWMEANCP.DF[,"bin"])) ){
      
      brmcp.mx <- matrix(data=NA, ncol=length(repeat.v), nrow=length(BINREP.MX[,"bins"]),
                         dimnames=list(NULL, repeat.v))
      brmcp.mx[, colnames(BINREP.MX[,-(1:3)])] <- BINREP.MX[,-(1:3)]
      
      brmcp.mx <- cbind(BINWMEANCP.DF[,c("wmeanCp0", "wmeanCp")], brmcp.mx)
      
      # Remove last bin because length <40kb
      brmcp.mx <- brmcp.mx[-length(brmcp.mx[,1]),]
      
      return(brmcp.mx)
      
    } else {
      stop(paste0("No data for ", chr, "."))
    }

  })
  
  return( do.call("rbind", chunk) )
  
}

elm.v <- read.csv(file=agerank.file, header=T, stringsAsFactors=F)[,"repName"]
elm.v.len <- length(elm.v)

toExport <- c("elm.v", "out.dir", "out.name")
P.LST <- foreach(itr=isplitVector(1:elm.v.len, chunks=nCPU), .inorder=T, .combine="c",
                 .export=toExport, .noexport=ls()[!ls()%in%toExport]
                           
) %op% {
  
  p.lst <- list()
  for(i in itr){
    
    elm <- elm.v[i]
    
    print(elm, quote=F) 
    
    pdf(file=paste0(out.dir, "/chrALL_", gcb, "_", elm, "_binrepcountsVsWmeanCp_scatter.pdf"), width=10, height=5)
    par(mfrow=c(1,2))
    
    for( wmCp in c("wmeanCp", "wmeanCp0") ){
      
      out.name <- paste0(wmCp, "_", elm)
      
      is.incl <- is.finite(BINREPMEANCP.MX[,wmCp]) & is.finite(BINREPMEANCP.MX[,elm])
      
      if(!inclbins0sites){
        # Consider only bins with at least 1 site, because bins without sites may be due to low copy number
        is.incl <- is.incl & BINREPMEANCP.MX[,elm] > 0
      } 
      
      incl.len <- sum(is.incl)
      
      plot.title <- paste0(out.name, "_chrAll_", gcb, "_", elm, "_", incl.len, "binsIncludedWith>0site")
      
      if(incl.len > 0){
        
        if(incl.len > 1){
          
          doCorTest(xval=BINREPMEANCP.MX[is.incl,wmCp], yval=BINREPMEANCP.MX[is.incl,elm],
                    alt="two.sided", exactpval=F, out.dir=out.dir, 
                    out.name=paste0(out.name, "_cortest.RData"))
          
          p.lst[[out.name]] <- try(
            makeHexbinggplot(xvar=BINREPMEANCP.MX[is.incl,wmCp], 
                             yvar=BINREPMEANCP.MX[is.incl,elm], 
                             bins=30, 
                             cuts=cuts,
                             xlab=wmCp,
                             ylab="element counts",
                             title=plot.title,
                             col=viridis(cuts))$hexplot
          )
          
        }
        
        plot(x=BINREPMEANCP.MX[is.incl,wmCp], y=BINREPMEANCP.MX[is.incl,elm], main=plot.title, cex.main=0.5,
             xlab=wmCp, ylab="element counts", col=adjustcolor("black", alpha.f=0.7))
        
      } else {
        print(paste0(out.name, ": No points to plot."))
        next
      }
    
    } # wmCp for loop end
    
    dev.off()
    
  } # itr for loop end
  
  return(p.lst)

}

p.arr <- ggarrange(plotlist=P.LST, nrow=5, ncol=5, legend=NULL)
ggexport(p.arr, height=50, width=50, 
         filename=paste0(out.dir, "/chrALL_", gcb, "_bins30_cuts", cuts, "_binrepcountsVsWmeanCp.pdf"))

# rm(list=ls()); gc()