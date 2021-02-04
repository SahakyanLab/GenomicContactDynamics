################################################################################
# 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/binmx/out_bindata_1perc_HiCNorm"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
meanCp.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binWeightedMeanCp")
mutCp.dir = paste0(wk.dir, "/out_mutCalcPerCp/WT_SEQ_rowSum")
out.dir = paste0(wk.dir, "/out_mutCalcVsMeanCp_bp/WT_SEQ_rowSum")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
src.id = "hg38ToHg19" # "Hg19" | "hg38ToHg19"
mut.v = c("All", "A>C", "A>G", "A>T", "C>A", "C>G", "C>T")
calc.v = c("Tmut", "Nmsite", "TmutDIVNmsite", "Nmsitenorm")
# Maximum number of bins with the same mean Cp0 and mean Cp value are 617 and 117
# bins only, respectively. Both are less than 1% of bins. So setting to frPerBoxplot
# any value greather than 0.01 is acceptable.
percPerBp = 5 # percPerBp of data in each boxplot
addjitter = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/splitNumericVector.R"))
source(paste0(lib, "/GG_bgr.R"))

makebp <- function(df, x="bin", y=calc, xlab="wmeanCp0-bin", ylab=calc, 
                   addjitter=addjitter, plot.id=paste0(mut.id, "_", src.id, "_", 
                                                       calc, "_", extCp)
                   ){
  
  binPerCp <- table(df[[x]])[levels(df[[x]])]
  binPerCp <- paste(paste(names(binPerCp), binPerCp, sep="="), collapse="_")
  
  # By default, ignore missing values in either the response or the group.
  eval(parse(text=paste0(
    'boxplot(', y, '~', x, ', outline=FALSE, data=df, xlab=xlab, 
     ylab=ylab, boxwex=0.6, cex.axis=1.2, col="#FDC776", cex.main=0.3,
     main=paste0(plot.id, "_", binPerCp))'
  )))
  
  if(addjitter){
    # Add data points
    levelprop.v <- summary(df[[x]])/nrow(df)
    for( i in 1:length(levels(df[[x]])) ){
      # Take the x-axis indices and add a jitter, proportional to the N in each level
      jitt <- jitter(rep(i, length(df[[y]])), amount=levelprop.v[i]/2)
      points(jitt, df[[y]], cex=1, col=adjustcolor("black",alpha.f=0.01), pch=16) 
      rm(jitt)
    }
  }
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
calc.v.len <- length(calc.v)*2
mut.v.len <- length(mut.v)
pdf(file=paste0(out.dir, "/", src.id, "_meanCp_percPerBp", percPerBp, "_boxplots.pdf"),
    height=mut.v.len*5, width=calc.v.len*5)
#png(file=paste0(out.dir, "/", src.id, "_meanCp_boxplots.png"),
#    height=mut.v.len*2*600, width=calc.v.len*2*600)
par(mfrow=c(mut.v.len, calc.v.len))

p.lst <- list()
for(mut in mut.v){
  
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=TRUE)
  
  # MUTCP.DF
  load(file=paste0(mutCp.dir, "/", mut.id, "_", src.id, "_mutCalcPerCp.RData"))
  rownames(MUTCP.DF) <- paste0(MUTCP.DF$chr, ".", MUTCP.DF$bins)
  MUTCP.DF$wmeanCp0 <- NA
  MUTCP.DF$wmeanCp <- NA
  chr.v <- unique(MUTCP.DF$chr)
  
  # Add mean Cp of bins in MUTCP.DF
  for(chr in chr.v){
    
    # BINWMEANCP.DF
    load(file=paste0(meanCp.dir, "/", gcb, "_", chr, "_weightedMeanCp.RData"))
    rownames(BINWMEANCP.DF) <- paste0(chr, ".", BINWMEANCP.DF$bin)
    
    if(!all(rownames(BINWMEANCP.DF)%in%rownames(MUTCP.DF))){
      stop(paste0(chr, ": Checkpoint 1."))
    }
    MUTCP.DF[rownames(BINWMEANCP.DF), "wmeanCp0"] <- BINWMEANCP.DF[rownames(BINWMEANCP.DF), "wmeanCp0"]
    MUTCP.DF[rownames(BINWMEANCP.DF), "wmeanCp"] <- BINWMEANCP.DF[rownames(BINWMEANCP.DF), "wmeanCp"]
    rm(BINWMEANCP.DF); gc()
    
  }

  wmeanCp0.TF <- is.finite(MUTCP.DF$wmeanCp0)
  wmeanCp.TF <- is.finite(MUTCP.DF$wmeanCp)
  
  for(calc in calc.v){
    
    #-------------------Weighted mean Cp (incl. Cp>=0 bins)
    bin.TF <- is.finite(MUTCP.DF[[calc]]) & wmeanCp0.TF
    df <- MUTCP.DF[bin.TF, c(calc, "wmeanCp0")]
    df <- df[order(df$wmeanCp0),]
    x <- stack( split( df[[calc]], ceiling(seq_along(df[[calc]])/(length(df[[calc]])*percPerBp/100)) ) )
    if( !identical(df[[calc]], x$values) ){ stop("Checkpoint 1.") }
    df <- cbind.data.frame(df, bin=x$ind)
    rm(x)
    
    extCp <- paste0("min=", min(df$wmeanCp0, na.rm=TRUE), "_",
                    "max=", max(df$wmeanCp0, na.rm=TRUE))
    
    makebp(df, xlab="wmeanCp0-bin", ylab=calc, addjitter=addjitter, 
           plot.id=paste0(mut.id, "_", src.id, "_", calc, "_", extCp))
    
    boundsSL <- c(min(df$wmeanCp0), ceiling(max(df$wmeanCp0)))
    slide <- splitNumericVector(x=df$wmeanCp0, d=1, action="SLIDE", numSLwind=20,
                                boundsSL=boundsSL)
    
    windowmid.v <- slide$midpoints; slide$midpoints <- NULL
    slide <- lapply(X=slide, FUN=function(ind){
      df[[calc]][ind]
    })
   
    df <- rbind(data.frame(stat="MEAN", stack(lapply(X=slide, FUN=mean))),
                data.frame(stat="MEDIAN", stack(lapply(X=slide, FUN=median)))
                )
    id <- paste0(mut.id, "_", src.id, "_", calc)
    p.lst[[id]] <- ggplot(data=df, aes(x=ind, y=values)) +
      geom_point(aes(colour=stat)) +
      labs(x="window", y=calc, 
           title=id) + 
      bgr2
    
    #-------------------Weighted mean Cp (incl. Cp>=1 bins)
    rm(df, bin.TF)
    
  }
  
  rm(MUTCP.DF, mut.id); gc()
  print(paste0(mut, " done!"), quote=FALSE)
  
} # mut.v for loop end

dev.off()

calc.v.len <- calc.v.len/2
p.arr <- ggarrange(plotlist=p.lst, nrow=mut.v.len, ncol=calc.v.len, legend=NULL)
ggexport(p.arr, height=mut.v.len*10, width=calc.v.len*10,
         filename=paste0(out.dir, "/", src.id, "_meanCp_percPerBp", percPerBp, "_scatter.pdf"))

# rm(list=ls()); gc()





