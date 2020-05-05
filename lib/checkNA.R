#check for missing values in datasets

################################################################################ 
checkNA <- function(x=dataset.mx){
  if( sum(apply(x, MARGIN=c(1,2), function(x) any(is.na(x))))!=0 ){
    stop("Missing values in datasets.")
  }
}
