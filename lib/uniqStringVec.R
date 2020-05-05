uniqStringVec <- function(vec=FEATURE.BIN.MX[FEATURE.BIN.MX$bin%in%bins.uniq, "name2"],
                          splitChar=";"){
  x <- unlist(x=strsplit(x=vec, split=splitChar), use.names=FALSE)
  return( paste(unique(x), collapse=splitChar) )
}