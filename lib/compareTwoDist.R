################################################################################
# Function to perform t test, Mann-Whitney test and Kolmogorov-Smirnoov test to 
# compare two distributions
### FUNCTION ###################################################################
compareTwoDist <- function(x, y){
  
  if( any(!is.finite(c(x,y))) ){
    
    warning("compareTwoDist(): Non-finite values in x and/or y.
            Removing those...")
    
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    
  }
  
  # Default var.equal = FALSE so no assumption of equal variances
  t <- t.test(x=x, y=y, alternative="two.sided", pool.sd=F, paired=F, conf.level=0.95)
  
  mw <- wilcox.test(x=x, y=y, alternative="two.sided", paired=F, conf.level=0.95)
  
  ks <- ks.test(x=x, y=y, alternative="two.sided")
  
  vals <- c(xmean=mean(x, na.rm=F), ymean=mean(y, na.rm=F),
            xmed=median(x, na.rm=F), ymed=median(y, na.rm=F),
            xmin=min(x, na.rm=F), ymin=min(y, na.rm=F),
            xmax=max(x, na.rm=F), ymax=max(y, na.rm=F)
            )
  vals <- format(x=vals, digits=5, trim=T, scientific=F)
  test.id <- paste0("\nmean_x;y=", vals[["xmean"]], ";", vals[["ymean"]],
                    "_median_x;y=", vals[["xmed"]], ";", vals[["ymed"]],
                    "_min_x;y=", vals[["xmin"]], ";", vals[["ymin"]],
                    "_max_x;y=", vals[["xmax"]], ";", vals[["ymax"]],
                    "\ntwosided_conflevel=0.95", 
                    "_ttest_", t$p.value, 
                    "_mwtest_", mw$p.value,
                    "_kstest_", ks$p.value)
  
  return( list(t=t, mw=mw, ks=ks, test.id=test.id) )
  
}

# rm(list=ls()); gc()