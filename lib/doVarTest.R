################################################################################
# Function to perform parametric one-way ANOVA and non-parametric counterpart
# Kruskal-Wallis H-test
### FUNCTION ###################################################################
doVarTest <- function(xval, grp, out.dir, out.name){
  
  if( any(is.na(c(xval, grp))) ){
    
    # Used stop to prompt user to manually check missing values that can mess up value-group
    # pairing
    stop("doANOVAKWH(): Non-finite values in xval and/or grp.")
    
  } else {
    
    ano <- aov(formula=xval~grp, data=data.frame(xval=xval, grp=grp))
    
    kwh <- kruskal.test(formula=xval~grp, data=data.frame(xval=xval, grp=grp))
    
    TEST <- list(ano=ano, kwh=kwh)
    save(TEST, file=paste0(out.dir, "/", out.name, "_varbasedtest.RData"))
    
  }
  
}

# rm(list=ls()); gc()