################################################################################
# Normalise contact value (i.e. complementarity values) based on distance between 
# contact bins (the farther the bins, that harder it should be to interact).
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(lib, "/simulation_lib/contactprobVsDistance.R"))
### FUNCTION ###################################################################
scaleContactByDist <- function(df.lst = 'list of two map values to be compared',
                               bin.len = 'resolution of bin',
                               out.filepath = 'output plot destination; no 
                               extension',
                               plot.title = 'plot title'){
  
  subj.TF <- grepl(x=names(df.lst), pattern="CII.cont", fixed=T) 
  
  if( sum(subj.TF)!=1 ){
    
    print("scaleContactByDist(): No map needing scaling.")
    return(df.lst)
    
  } else {
    
    subj <- names(df.lst)[subj.TF]
    ref <- names(df.lst)[!subj.TF]
    
  }
  
  # Get average contact value per gap in bins (i.e. per Hi-C diagonal)
  # Note that contactprobVsDistance() includes 0 values
  
  refgap.meanval.df <- contactprobVsDistance(df=df.lst[[ref]],
                                             bin.len=bin.len, scale.diag.ind=NA,
                                             smooth.method="NONE")$bydiagonal[,c("diag.ind", "val.ave")] 
  
  colnames(refgap.meanval.df) <- c("gapjMINUSi", "mean.val.nonNAs")
  refgap.meanval.df$gapjMINUSi <- as.numeric(as.character(refgap.meanval.df$gapjMINUSi))
  rownames(refgap.meanval.df) <- refgap.meanval.df$gapjMINUSi
  
  # Normalise
  max.val <- max(refgap.meanval.df$mean.val.nonNAs, na.rm=F)
  min.val <- min(refgap.meanval.df$mean.val.nonNAs, na.rm=F)
  refgap.meanval.df$mean.val.nonNAs.norm <- (refgap.meanval.df$mean.val.nonNAs - min.val) / (max.val-min.val)
  #refgap.meanval.df$mean.val.nonNAs.norm <- refgap.meanval.df$mean.val.nonNAs / max.val
  
  # Fit non-linear function to raw ref values
  #x.val <- refgap.meanval.df$gapjMINUSi 
  #y.val <- refgap.meanval.df$mean.val.nonNAs
  #ref.fit <- nls(formula=y.val ~ a*(x.val^b))
  #refgap.meanval.df$pred <- ref.fit$m$predict()
  
  ## Fit non-linear function to normalised ref values
  x.val <- refgap.meanval.df$gapjMINUSi 
  y.val <- refgap.meanval.df$mean.val.nonNAs.norm
  ref.fit <- nls(formula=y.val ~ a*(x.val^b))
  refgap.meanval.df$pred <- ref.fit$m$predict()
  
  ## Normalise pred values
  #max.val <- max(refgap.meanval.df$pred, na.rm=F)
  #min.val <- min(refgap.meanval.df$pred, na.rm=F)
  #refgap.meanval.df$pred <- (refgap.meanval.df$pred - min.val) / (max.val-min.val)
  ##refgap.meanval.df$pred <- refgap.meanval.df$pred / max.val
  
  rm(x.val, y.val)
  
  # Transform CII continuous to positive values that can be transformed
  # to probabilities. Did eumiro's answer so distribution will only be translated:
  # https://stackoverflow.com/questions/3931419/turn-a-negative-number-into-a-positive-for-probability
  
  nonNA.TF <- !is.na(df.lst[[subj]]$value)

  # +1 so probability for minimum value will be 1 (not 0)
  df.lst[[subj]]$value[nonNA.TF] <- df.lst[[subj]]$value[nonNA.TF] - 
                                    min(df.lst[[subj]]$value, na.rm=T) + 1
  
  # Multiply subj values by ref pred per gap given by fitted function
  
  subj.gapbin.v <- abs( df.lst[[subj]]$j-df.lst[[subj]]$i )
  
  is.subjgapwithMult <- subj.gapbin.v%in%refgap.meanval.df$gapjMINUSi
  subj.mult.v <- rep(x=NA, times=length(subj.gapbin.v))
  subj.mult.v[is.subjgapwithMult] <- refgap.meanval.df[ as.character(subj.gapbin.v[is.subjgapwithMult]),
                                                        "pred" ]
  
  df.lst[[subj]]$value <- df.lst[[subj]]$value * subj.mult.v
  
  #
  subjgap.meanval.df <- contactprobVsDistance(df=df.lst[[subj]],
                                              bin.len=bin.len, scale.diag.ind=NA,
                                              smooth.method="NONE")$bydiagonal[,c("diag.ind", "val.ave")] 
  
  colnames(subjgap.meanval.df) <- c("gapjMINUSi", "mean.val.nonNAs")
  subjgap.meanval.df$gapjMINUSi <- as.numeric(as.character(subjgap.meanval.df$gapjMINUSi))
  rownames(subjgap.meanval.df) <- subjgap.meanval.df$gapjMINUSi
  
  max.val <- max(subjgap.meanval.df$mean.val.nonNAs, na.rm=F)
  min.val <- min(subjgap.meanval.df$mean.val.nonNAs, na.rm=F)
  subjgap.meanval.df$mean.val.nonNAs.norm <- (subjgap.meanval.df$mean.val.nonNAs - min.val) / (max.val-min.val)
  #subjgap.meanval.df$mean.val.nonNAs.norm <- subjgap.meanval.df$mean.val.nonNAs / max.val

  # 
  png(file=paste0(out.filepath, ".png"), width=300*20, height=300*10, res=300)
  par(mfrow=c(1,2))
 
  plot(density(df.lst[[subj]]$value, na.rm=T), 
       main=paste0("subj", subj, "values_multplied by ref", ref, " mean.val.nonNAs divided by max value"),
       cex.main=0.5, col=adjustcolor(col="darkred"), xlim=c(0,4))
  
  plot.id <- paste0("refinblue", ref, "_subjinred", subj, 
                    "_raw=brightshade_norm=darkshade_fittedfunction=black", 
                    "_fittedparameters_",
                    paste(paste(names(ref.fit$m$getPars()), ref.fit$m$getPars(), sep="="), collapse="_"),
                    "\nrefrawVSsubjrawscaled_cor_pearson=",
                    cor(x=refgap.meanval.df$mean.val.nonNAs, 
                        y=subjgap.meanval.df$mean.val.nonNAs, method="pearson"),
                    "_spearman=",
                    cor(x=refgap.meanval.df$mean.val.nonNAs, 
                        y=subjgap.meanval.df$mean.val.nonNAs, method="spearman"),
                    "_gapResolution=", bin.len, "bp_",
                    "refmingap=", min(refgap.meanval.df$gapjMINUSi),
                    "refmaxgap=", max(refgap.meanval.df$gapjMINUSi),
                    "subjmingap=", min(subjgap.meanval.df$gapjMINUSi),
                    "subjmaxgap=", max(subjgap.meanval.df$gapjMINUSi)
                    )
  
  plot(x=refgap.meanval.df$gapjMINUSi, 
       y=refgap.meanval.df$mean.val.nonNAs, 
       cex=0.5, cex.main=0.5, col=adjustcolor(col="blue", alpha.f=0.5), 
       xlab="gapjMINUSi", ylab="mean.val.nonNAs divided by max value",
       main=paste0(plot.title, "\n", plot.id, "_scaleContactByDistPlot"), ylim=c(0,2))
  
  points(x=refgap.meanval.df$gapjMINUSi, 
         y=refgap.meanval.df$mean.val.nonNAs.norm, 
         cex=0.5, col=adjustcolor(col="darkblue", alpha.f=0.5))
  
  points(x=refgap.meanval.df$gapjMINUSi, 
         y=refgap.meanval.df$pred, 
         cex=0.5, col=adjustcolor(col="black", alpha.f=0.5))
  
  points(x=subjgap.meanval.df$gapjMINUSi, 
         y=subjgap.meanval.df$mean.val.nonNAs,#- (max(subjgap.meanval.df$mean.val.nonNAs, na.rm=T) - 1), 
         cex=0.5, col=adjustcolor(col="red", alpha.f=0.5))
  
  points(x=subjgap.meanval.df$gapjMINUSi, 
         y=subjgap.meanval.df$mean.val.nonNAs.norm, 
         cex=0.5, col=adjustcolor(col="darkred", alpha.f=0.5))

  dev.off()
  
  return(df.lst)
  
}

# rm(list=ls()); gc()