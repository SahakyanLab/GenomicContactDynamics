############################################################################### 
# Calculate center of mass
###############################################################################
getCOM <- function(weight = weight.v,
                   coordinates = coord.df){
  
  # Make sure coordinates and weights were converted to numeric 
  # in the getcmmXYZradi function
  Totweight <- sum(weight)
  COM <- unlist(lapply(X=coordinates, FUN=function(coord){
                sum(coord*weight)/Totweight
              }), use.names=TRUE
             )
  return(COM)
  
}
  