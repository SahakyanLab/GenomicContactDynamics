getRadDist <- function(xyz1, xyz2) {
  
  dist <- sqrt( ( (xyz1[1] - xyz2[1])^2 ) +
                ( (xyz1[2] - xyz2[2])^2 ) +
                ( (xyz1[3] - xyz2[3])^2 )   )
  return(dist)
  
}