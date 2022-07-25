################################################################################
# Calculate consensus value for contact out of values from contacting bins
### FUNCTION ###################################################################
ijconsensusFUN <- function(i, j, METHOD){
  
  M.v <- tolower(unlist(strsplit(x=METHOD, split=".", fixed=T)))
  
  if(M.v[1] != "c"){ # Stat per bin
  
    eval(parse(text=paste0(
      "i <- ", M.v[1], "(i, na.rm=T)"
    )))
    
    eval(parse(text=paste0(
      "j <- ", M.v[1], "(j, na.rm=T)"
    )))
    
  } # else, pool i and j values
  
  eval(parse(text=paste0(
    "consensus <- ", M.v[2], "(c(i,j), na.rm=T)"
  )))
  
  return(consensus)
  
}
################################################################################

# rm(list=ls()); gc()