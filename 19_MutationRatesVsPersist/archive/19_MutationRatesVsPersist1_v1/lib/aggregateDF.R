############################################################################### 
# Specific function to aggregate boxplot-style dataframe using mean and median
###############################################################################
aggregateDF <- function(ind=x$Cp, values=x[[calc]]
                        ){
  
  mean.df <- aggregate(x=values, by=list(ind), FUN=mean, na.rm=TRUE)
  mean.df$Group.1 <- as.numeric(mean.df$Group.1)
  mean.df$id <- "MEAN"
  rownames(mean.df) <- paste0(mean.df$Group.1, "_MEAN")
  
  med.df <- aggregate(x=values, by=list(ind), FUN=median, na.rm=TRUE)
  med.df$Group.1 <- as.numeric(med.df$Group.1)
  med.df$id <- "MEDIAN"
  rownames(med.df) <- paste0(med.df$Group.1, "_MEDIAN")
  
  df <- rbind.data.frame(mean.df, med.df)
  colnames(df) <- c("ind", "values", "id")
  
  return(df)
  
}
