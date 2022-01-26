## FUNCTION ####################################################################
# A function to extract persistent chromatine contacts from supplied MELT.MX
# object (its upper.tri component).
################################################################################
getMELTMXpersist <- function(MELT.MX   = MELT.MX,
                            # Minimum allowed genomic distance between contacts:
                             min.dist  = 2000000,  # 2Mbp
                         # The resolution of Hi-C contact matrix behind MELT.MX:
                             hic.resol = 40000,    # 40kbp
              # The percentage showing that the value cutoff should be below the
              # value of the non-0 and >=min.dist contacts in each individual
              # dataset (tissue/cell type):  
                             top.perc  = 100,          # 100% will consider all.
                        # Minimum number of tissue/cell types (out of 21), where
                        # the top.perc should be kept (21 for all).
                             min.tiss  = 1
                             #control   = TRUE,
                             #control.interval.bin = 10,
                             #seed = 845
                            ){
  
  MELT.MX$upper.tri <- MELT.MX$upper.tri[((MELT.MX$upper.tri[,"j"]-
                                           MELT.MX$upper.tri[,"i"]-1) >=
                                           (min.dist/hic.resol)),]
  
  MELT.MX$upper.tri.nocontact <-
           MELT.MX$upper.tri.nocontact[((MELT.MX$upper.tri.nocontact[,"j"]-
                                         MELT.MX$upper.tri.nocontact[,"i"]-1) >=
                                         (min.dist/hic.resol)),]
  
  TrueFalse.MELT.MX <- MELT.MX$upper.tri
  
  tis.len <- length(MELT.MX$upper.tri[1,])
  for(tis in 3:tis.len){
  #for(tis in 3:23){
    value <- MELT.MX$upper.tri[,tis]
    value <- sort(value[value!=0], decreasing=TRUE)
    value <- value[floor(length(value)*top.perc/100)] # Should be higher or
    # equal than <value> for satisfying the top.perc cutoff.
    TrueFalse.MELT.MX[,tis] <- (MELT.MX$upper.tri[,tis]>=value)
  }
  
  #rowsums.TrueFalse.MELT.MX <- rowSums(TrueFalse.MELT.MX[,3:23])
  rowsums.TrueFalse.MELT.MX <- rowSums(TrueFalse.MELT.MX[,3:tis.len, drop=FALSE])
  #rowsums.MELT.MX.uppertri  <- rowSums(MELT.MX$upper.tri[,3:23])
  rowsums.MELT.MX.uppertri  <- rowSums(MELT.MX$upper.tri[,3:tis.len, drop=FALSE])
  hits.ids <- which( rowsums.TrueFalse.MELT.MX >= min.tiss )
    
  length.left <- length(hits.ids)
  if(length.left!=0){
    print(paste0("getMELTMXpersist: ",length.left,
            " entries remained after the provided value cuttof that occured in at least ",
                 top.perc,"% data, in at least ",min.tiss," tissues."),quote=FALSE)
  } else {
    stop(paste0("getMELTMXpersist: no entry remained from the provided ",
    top.perc,"% cutoff in at least ",min.tiss," tissues."))
  }

  RESULT <- NULL
  RESULT$hits    <- MELT.MX$upper.tri[hits.ids,]
  # tissues num. the conditions are true
  RESULT$ntis    <- as.vector(rowsums.TrueFalse.MELT.MX[hits.ids])
  # sum of all contact values across all tissues
  RESULT$valsum  <- as.vector(rowsums.MELT.MX.uppertri[hits.ids]) 
  RESULT$control <- MELT.MX$upper.tri.nocontact
  
    # ToDo: matching control subsampled distribution structure to the source.
    #if(control){
    # # Contact intervals from the hits (as names), and their count (as values).
    # fltrd.cntct.interv <- table(RESULT$hits[,"j"]-RESULT$hits[,"i"])
    # #  51  52  53 ... # names
    # # 33  35  28  ... # number of occurrence
    # fltrd.cntct.interv.names <- as.numeric(names(fltrd.cntct.interv))
    # interval.range   <- range(fltrd.cntct.interv.names)
    # interval.seq     <- seq(from=interval.range[1],
    #                         to=interval.range[2], by=control.interval.bin)
    # interval.seq.len <- length(interval.seq)
    # # 51   61   71   81   91  101  111  121 ...  
    # 
    # nocontact.intervals <- (MELT.MX$upper.tri.nocontact[,"j"]-
    #                         MELT.MX$upper.tri.nocontact[,"i"])
    # 
    # for(fci in interval.seq){
    #   interval <- c(fci, fci+control.interval.bin-1)
    #   
    #   ind1 <- which( nocontact.intervals >= interval[1] &
    #                  nocontact.intervals  < interval[2] )
    #   
    #   how.many <- sum(
    #       fltrd.cntct.interv[which(fltrd.cntct.interv.names >= interval[1] &
    #                                fltrd.cntct.interv.names  < interval[2])]
    #                  )
    #     x[!(x %in% y)]
    # }
    #}
  return(RESULT)  
}
################################################################################









