############################################################################### 
# Specific function to aggregate boxplot-style dataframe using mean and median
###############################################################################
aggregateDF <- function(ind, values, FUN){
  
  if( !is.factor(ind) ){
    stop("aggregateDF: ind not a factor as expected.")
  }
  # aggregate() coerces by as factor
  eval(parse(text=paste0(
    'df <- aggregate(x=values, by=list(ind), FUN=', FUN, ', na.rm=T)'
  )))
  
  df$Group.1 <- as.character(df$Group.1)
  colnames(df) <- c("ind", "values")
  
  return(df)
  
}
