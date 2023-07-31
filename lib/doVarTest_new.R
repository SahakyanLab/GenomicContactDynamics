################################################################################
# Function to perform parametric one-way ANOVA and non-parametric counterpart
# Kruskal-Wallis H-test
# Add two-way ANOVA for unbalanced designs
# Refer to http://www.sthda.com/english/wiki/two-way-anova-test-in-r#compute-two-way-anova-test-in-r-for-unbalanced-designs
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#library(car) # Anova(), ANOVA for unbalanced dataset
### FUNCTION ###################################################################
doVarTest <- function(xval, grp, out.dir, out.name, lightenAOV = F){
  
  if( any(is.na(c(xval, grp))) ){
    
    # Used stop to prompt user to manually check missing values that can mess up value-group
    # pairing
    stop("doVarTest(): Non-finite values in xval and/or grp.")
    
  } else {
    
    ano <- aov(formula=xval~grp, data=data.frame(xval=xval, grp=grp))
  
    #ano_unbTyp3 <- NULL
    #if(do2WAYANOVA){ ano_unbTyp3 <- Anova(ano, type="III") } 
    
    if(lightenAOV){ ano <- summary(ano) }
    
    #
    kwh <- kruskal.test(formula=xval~grp, data=data.frame(xval=xval, grp=grp))
    
    #
    TEST <- list(ano=ano, kwh=kwh) #, ano_unbTyp3=ano_unbTyp3) # kwh[c("statistic", "p.value")]
    
    if( !is.null(out.dir) & !is.null(out.name) ){
      save(TEST, file=paste0(out.dir, "/", out.name, "_varbasedtest.RData"))
    } 
    
    return(TEST)
    
  }
  
}

# rm(list=ls()); gc()