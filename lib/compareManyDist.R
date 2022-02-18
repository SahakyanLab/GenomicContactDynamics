################################################################################
# Function to perform pairwise comparison of multiple distributions using t-test 
# & Mann-Whitney test.
### FUNCTION ###################################################################
compareManyDist <- function(xval, grp, #'grp should be a factor' 
                            alt="two.sided", out.dir, out.name){
  
  if( any(is.na(c(xval, grp))) ){
    
    # Used stop to prompt user to manually check missing values that can mess up value-group
    # pairing
    stop("compareManyDist(): Non-finite values in xval and/or grp.")
    
  } else {
    
    # pool.sd = F -> no assumption of equal variances
    pt <- pairwise.t.test(x=xval, g=grp, p.adjust.method="BH", alternative=alt, 
                          pool.sd=F, paired=F, conf.level=0.95)
    
    pmw <- pairwise.wilcox.test(x=xval, g=grp, p.adjust.method="BH", alternative=alt,
                                paired=F, conf.level=0.95)
    
    meanmed <- by( data=xval, INDICES=list(grp), 
                   FUN=function(x) c(MEAN=mean(x), MED=median(x)) )
   
    TEST <- list(pt=pt, pmw=pmw, meanmed=meanmed, alt=alt)
    save(TEST, file=paste0(out.dir, "/", out.name, "_pairwisedifftest.RData"))
    
  }
  
}

# rm(list=ls()); gc()