################################################################################
findHPM <- function(x=x){
  x <- gsub("\\-|\\[|\\]|\\d", "", x, ignore.case=TRUE)
  x <- strsplit(x, split=";")[[1]]
  total <- length(x)
  mammal.ind <- grep("Cjac|Ocun|Mmus|Rnor|Ecab|Btau|Sscr|Cfam|Oari|Fcat", x, 
                     ignore.case=TRUE)
  mammal.count <- length(mammal.ind)
  x <- x[-mammal.ind]
  remove.hsap <- gsub("h|s|a|p", "", 
                    x, ignore.case=TRUE)
  human.count <- length(which((nchar(remove.hsap))==0))
  primate.count <- total-(human.count+mammal.count)
  if(total!=(human.count+primate.count+mammal.count)){
    stop("HMP counts not adding up.")
  }
  f <- as.numeric( c(human.count/total, primate.count/total, mammal.count/total, total) )
  names(f) <- c("fH", "fP", "fM", "total")
  f
}