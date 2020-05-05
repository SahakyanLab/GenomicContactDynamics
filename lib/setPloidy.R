setPloidy <- function(domainObj = DOMAIN.df, ploidy = 2,
                      chrCol = "chr", idCol = "id"){

  n.row <- nrow(domainObj)
  #make bead id according to ploidy
  lab.v <- paste(chr <- rep(x=domainObj[[chrCol]], each=ploidy), 
                 rbind(rep("A", n.row), rep("B", n.row)), sep="_")
  coord.v <- as.vector(do.call("rbind", strsplit(domainObj[[idCol]], split=":"))[,2])
  id <- paste(lab.v, rep(coord.v, each=ploidy), sep=":")
  
  #duplicate DomainVsFeature 
  domainObj <- domainObj[rep(x=row.names(domainObj), each=ploidy), 
                         1:ncol(domainObj)]
  domainObj[[idCol]] <- id
  
  return(domainObj)
}
