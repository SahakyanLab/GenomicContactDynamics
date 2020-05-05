makeScinote <- function(x){
  if(grepl(pattern="e", x=x)){
    return( gsub(pattern="e", x=x, replacement="%*%10^") )
  } else {
    return(x)
  }
}