#RecombinationRatesVsPersist

#REMINDERS BEFORE RUNNING
#don't forget to modify the included rows of PERSIST.MX, adjusted for testing only 

#recombination rate file structure: chr-position(bp)-rate(cM/Mb)-map(cM)

################################################################################
#Directories

#functions to source
#lib <- "/Users/ltamon/ProjectBoard/lib"
lib <- "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/2_RecombinationRatesVsPersist"
objective.dir <- "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/2_RecombinationRatesVsPersist"
#recombination rates data (Myers et al., 2005; Munch et al., 2014)
#recomRates.dir <- "/Users/ltamon/Database/recomRates_2011-01_phaseII_B37_Myers"
recomRates.dir <- "/t1-data/user/ltamon/Database/recomRates_2011-01_phaseII_B37_Myers"
#HiC_Human21 persistent contact files
#persist.dir <- "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
persist.dir <- "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
#dir.create(paste0(objective.dir, "/out_RecomVsPersist"))
output.dir <- paste0(objective.dir, "/out_RecomVsPersist") 

################################################################################
#Packages

##CRAN
library("foreach")
library("data.table") 

##BIOCONDUCTOR
library("IRanges") 

################################################################################
#Set values 

#length of ij bins
bin.length <- 40000L

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gc.v <- c(2, "05") 

#no recombination data for Y and MT from Myers2005/Munsch2014
chr.v <- c(1:22, "X")

#column names of contacting bins in the persistent contact matrices
i <- "i"
j <- "j"

################################################################################ 
#Functions

#source(paste0(lib, "/UTL_NCPU.R"))
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))

################################################################################
################################################################################
#PART1 - Assign recombination rates to bins

#nested foreach to get data for each chromosome for each gap 
foreach(gc=gc.v, .inorder=TRUE) %do% {
  foreach(chr=chr.v, .inorder=TRUE) %do% {
    recomRates.mx <- fread(file=paste0(recomRates.dir, "/genetic_map_GRCh37_chr", chr,".txt"), 
                           header=TRUE, data.table=FALSE)
    colnames(recomRates.mx) <- c("Chr", "POSbp", "RATEcMperMb", "MAPcM")
    
    #if( sum(apply(recomRates.mx, MARGIN=c(1,2), function(x) any(is.na(x))))!=0 ){
    #  stop("Missing values in recombination rate file.")}
  #} foreach closed to check if recombination rate files have missing values,
    #all files - no NA
  
    load(paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))
    persist.mx   <- PERSIST.MX$hits[,c(i, j)]
    rownames(persist.mx) <- NULL
    persist.ntis <- PERSIST.MX$ntis   
   
    #get the relevant bins (from persist matrix) for that chromosome
    #bins sorted to increasing values so the ranges are also increasing
    bins.uniq <- sort(unique(unlist(persist.mx[, c(i, j)])), decreasing=FALSE)
    bin.end   <- bins.uniq*bin.length
    bin.start <- bin.end-bin.length+1
    
    #assign positions of recombination rates to bins of the chromosome
    #query   <- recombination rates (to be classified)
    #subject <- bins (classification)
    olap.mx <- WhichOverlap(start.query=recomRates.mx[,"POSbp"], 
                            end.query=recomRates.mx[,"POSbp"], 
                            space.query=rep("a", nrow(recomRates.mx)),
                            start.subject=bin.start, 
                            end.subject=bin.end, 
                            space.subject=rep("a",length(bin.end)),
                            maxgap=0L, minoverlap=1L) 
    
    #make the matrix with the recombination rates assigned to bins
    #for bins with multiple rates, get mean, median and sd
    rates.bin      <- cbind(binAssignment=bins.uniq[olap.mx[,"subject"]], 
                            recomRates.mx[olap.mx[,"query"],])
    count.indv.r.b <- by(rates.bin[,"RATEcMperMb"], rates.bin[,"binAssignment"], 
                         function(x) length(x))
    mean.r.b       <- by(rates.bin[,"RATEcMperMb"], rates.bin[,"binAssignment"], 
                         function(x) mean(x))
    median.r.b     <- by(rates.bin[,"RATEcMperMb"], rates.bin[,"binAssignment"], 
                         function(x) median(x))
    #sd can give NA values when n=1
    sd.r.b         <- by(rates.bin[,"RATEcMperMb"], rates.bin[,"binAssignment"], 
                         function(x) sd(x))
    indv.pos.b     <- by(rates.bin[,"POSbp"], rates.bin[,"binAssignment"], 
                         function(x) paste(x, collapse = ","))
    indv.r.b       <- by(rates.bin[,"RATEcMperMb"], rates.bin[,"binAssignment"], 
                         function(x) paste(x, collapse = ","))
    indv.map.b     <- by(rates.bin[,"MAPcM"], rates.bin[,"binAssignment"], 
                         function(x) paste(x, collapse = ","))
    
    bin.assigned <- unique(rates.bin[,"binAssignment"])
    
    RECOM.BIN.MX <- cbind(Bin=bin.assigned, 
                          MEANrb=mean.r.b, 
                          MEDIANrb=median.r.b, 
                          SDEVrb=sd.r.b, 
                          COUNTINDVrb=count.indv.r.b,
                          INDVRATEScMperMb=indv.r.b, 
                          INDVPOSbp=indv.pos.b, 
                          INDVMAPcM=indv.map.b) 
    row.names(RECOM.BIN.MX) <- NULL
    
#---------------------------------------
#PART2 - Make a big matrix with the ij contacts, recombination rates 
#and relevant statistical values
   
    #isolate only contacts that have recombination rate 
    #values for both i and j
    ind.and  <- which(persist.mx[,i]%in%bin.assigned & persist.mx[,j]%in%bin.assigned)
    ntis.and <- persist.ntis[ind.and] 
    
    #subset of PERSIST.MX$hits with recombination rates for i AND j bins
    i.and          <- persist.mx[ind.and,i]
    j.and          <- persist.mx[ind.and,j]
    
    #assign the subset of i and j bins to recombination rates from RECOM.BIN.MX
    #query   <- i or j bins with recombination rates
    #subject <- bin from RECOM.BIN.MX (stored as bin.assigned variable), 
    #which has recombination data (etc.) for a particular bin
    olap.i.and <- WhichOverlap(start.query=i.and, 
                               end.query=i.and, 
                               space.query=rep("a", length(i.and)),
                               start.subject=bin.assigned, 
                               end.subject=bin.assigned, 
                               space.subject=rep("a", length(bin.assigned)),
                               maxgap=0L, minoverlap=1L)
    olap.j.and <- WhichOverlap(start.query=j.and, 
                               end.query=j.and, 
                               space.query=rep("a", length(j.and)),
                               start.subject=bin.assigned, 
                               end.subject=bin.assigned, 
                               space.subject=rep("a", length(bin.assigned)),
                               maxgap=0L, minoverlap=1L)
    
    ij.and.mx <- cbind(RECOM.BIN.MX[olap.i.and[,"subject"], 
                                    c("MEANrb", "MEDIANrb", "INDVRATEScMperMb")], 
                       RECOM.BIN.MX[olap.j.and[,"subject"], 
                                    c("MEANrb", "MEDIANrb", "INDVRATEScMperMb")])
                       
    ind.mean      <- which(colnames(ij.and.mx)=="MEANrb")
    ind.med       <- which(colnames(ij.and.mx)=="MEDIANrb")
    ind.indvrates <- which(colnames(ij.and.mx)=="INDVRATEScMperMb")
   
    #statistical analyses pooling individual rates of i and j bins
    result.pool  <- lapply(list(INDVRATEScMperMb=ij.and.mx[,ind.indvrates]),
      function(x){
        ij.pool        <- apply(x, MARGIN=1, function(x) paste(x, collapse=","))
        ij.pool        <- lapply(ij.pool, function(x) unlist(strsplit(x, split=",")) )
        ij.pool        <- lapply(ij.pool, function(x) as.numeric(x))
        mean.ij.pool   <- lapply(ij.pool, function(x) mean(x))
        median.ij.pool <- lapply(ij.pool, function(x) median(x))
        sd.ij.pool     <- lapply(ij.pool, function(x) sd(x)) #NA when n=1
        min.ij.pool    <- lapply(ij.pool, function(x) min(x))
        max.ij.pool    <- lapply(ij.pool, function(x) max(x))
        all.val <- cbind(i=as.character(i.and), j=as.character(j.and),
                         MEANijP= mean.ij.pool , MEDIANijP=median.ij.pool, 
                         SDEVijP=sd.ij.pool, MINijP=min.ij.pool, MAXijP=max.ij.pool)
      })
      
    #statistical analyses for the values of contacting i and j bins
    #mean, median, sd, min, max, difference
    mx.list <- list()
    mx.list <- list(MEANrb=ij.and.mx[,ind.mean], MEDIANrb=ij.and.mx[,ind.med])
    result  <- lapply(mx.list, function(x){
      x         <- apply(x, MARGIN=2, as.numeric)
      mean.ij   <- apply(x, MARGIN=1, mean)
      median.ij <- apply(x, MARGIN=1, median)
      sd.ij     <- apply(x, MARGIN=1, sd) #NA when n=1
      min.ij    <- apply(x, MARGIN=1, min)
      max.ij    <- apply(x, MARGIN=1, max)
      diff.ij   <- abs(x[,1]-x[,2])
      all.val <- cbind(i=as.character(i.and), j=as.character(j.and), 
                       MEANij=mean.ij, MEDIANij=median.ij, 
                       SDEVij=sd.ij, MINij=min.ij, MAXij=max.ij, DIFFij=diff.ij)
    })
    
    incl  <- c()
    total <- c()
    perc  <- c()
    rates.stat <- c(incl=nrow(olap.mx), tot=nrow(recomRates.mx), 
                    perc=(nrow(olap.mx)/nrow(recomRates.mx))*100)
    ij.stat    <- c(incl=nrow(ij.and.mx), tot=nrow(persist.mx), 
                    perc=(nrow(ij.and.mx)/nrow(persist.mx)*100))
    
    RECOM.PERSIST <- list(BINVAL=result,
                          POOL=result.pool$INDVRATEScMperMb,
                          NTIS=ntis.and,
                          RECOM.BIN.MX=cbind(Chr=rep(paste0("chr", chr), length(indv.pos.b)),
                                             RECOM.BIN.MX),
                          RATESstat=round(rates.stat, digits=2), 
                          CONTACTSstat=round(ij.stat, digits=2)
    )
    
    save(RECOM.PERSIST, file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_RecomVsPersist.RData"))
    
  } 
} 

#rm(list=ls())

