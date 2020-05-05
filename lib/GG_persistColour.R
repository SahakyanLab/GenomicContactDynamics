################################################################################
PersistScoreColour <- function(PersistScores = ntisUsed){
  hexColour <- lapply(PersistScores, function(x) {
    rgb.v <- c(255,255,255)-(12*x)
    rgb(rgb.v[1], rgb.v[2], 255, max=255)
  } )
   unlist(hexColour)
}
