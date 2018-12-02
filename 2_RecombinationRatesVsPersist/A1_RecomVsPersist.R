#RecombinationRatesVsPersist

#REMINDERS BEFORE RUNNING
#don't forget to modify the included rows of PERSIST.MX, adjusted for testing only 

#recombination rate file structure: chr-position(bp)-rate(cM/Mb)-map(cM)

#for testing, make sure to include at least 500 of ij contacts from the persist
#matrix to make sure that you include ij contacts both having recombination 
#rates values
#i still have to find a way to do NEXT (skip iteration) in foreach 
#to account for this
#if(length(ind.and)==0) {
#  break # skip current iteration and go to next iteration
#  noOUTPUT[paste0("chr", chr)] <- gc
#} 

################################################################################
#Directories

#functions to source
lib <- "/Users/ltamon/ProjectBoard/lib"
#main directory for the task
objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/2_RecombinationRatesVsPersist"
#recombination rates data (Myers et al., 2005; Munch et al., 2014)
recomRates.dir <- "/Users/ltamon/Database/recomRates_2011-01_phaseII_B37_Myers"
#HiC_Human21 persistent contact files
persist.dir <- "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
#dir.create(paste0(objective.dir, "/RecomVsPersist"))
output.dir <- paste0(objective.dir, "/RecomVsPersist")  

################################################################################
#Packages

##CRAN
#install.packages("foreach") 
library("foreach")
#install.packages("data.table") 
library("data.table") #for fread

##BIOCONDUCTOR
#source("https://bioconductor.org/biocLite.R")
#biocLite("IRanges")
library("IRanges") #for WhichOverlap

################################################################################
#Set values 

NCPU <- 1

#length of ij bins
bin.length <- 40000L

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gc.v <- c(2, "05") 

#no recombination data for Y and MT from Myers2005/Munsch2014
chr.v <- c(3,4) #c(1:22, "X")

#column names of contacting bins in the persistent contact matrices
i <- "i"
j <- "j"

################################################################################ 
#Functions

source(paste0(lib, "/UTL_NCPU.R"))
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))

################################################################################
################################################################################
#PART1 - Assign recombination rates to bins

#nested foreach to get data for each chromosome for each gap 
foreach(gc=gc.v, .inorder=FALSE) %op% {
  
  allMEANrb     <- list()
  allMEDIANrb   <- list()
  allrates.stat <- list()
  allij.stat    <- list()
 
  foreach(chr=chr.v, .inorder=FALSE) %op% {
    #read recombination rate file 
    recomRates.mx <- fread(file=paste0(recomRates.dir, "/genetic_map_GRCh37_chr", chr,".txt"), 
                           header=TRUE, data.table=FALSE)
    colnames(recomRates.mx) <- c("Chr", "POSbp", "RATEcMperMb", "MAPcM")
    
    #if( sum(apply(recomRates.mx, MARGIN=c(1,2), function(x) any(is.na(x))))!=0 ){
    #  stop("Missing values in recombination rate file.")}
  #} foreach closed to check if recombination rate files have missing values, all files - no NA
  
    #load persistent contact files 
    load(paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))
    #PERSIST.MX$hits and PERSIST.MX$ntis are in the same order
    persist.mx   <- PERSIST.MX$hits[1:5000,]  
    persist.ntis <- PERSIST.MX$ntis[1:5000]   
   
    #get the relevant bins (from persist matrix) for that chromosome
    #bins sorted to increasing values so the ranges are also increasing
    bins.uniq <- sort(unique(unlist(persist.mx[, c(i, j)])), decreasing=FALSE)
    #end points of 40-Kb bins
    bin.end   <- bins.uniq*bin.length
    #start points of 40-Kb bins
    bin.start <- bin.end-bin.length+1
    
    #assign positions of recombination rates to bins of the chromosome
    #findoverlaps based on actual ranges of bins
    #Note that WhichOverlap(IRanges) returns indices
    #query   <- recombination rates (to be classified)
    #subject <- bins (classification)
    olap.mx <- WhichOverlap(start.query=recomRates.mx[,"POSbp"], 
                            end.query=recomRates.mx[,"POSbp"], 
                            space.query=rep("a", nrow(recomRates.mx)),
                            start.subject=bin.start, 
                            end.subject=bin.end, 
                            space.subject=rep("a",length(bin.end)),
                            maxgap=-1L, minoverlap=1L) #default
    
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
    
    RECOM.BIN.MX <- cbind(Chr=rep(paste0("chr", chr), length(indv.pos.b)), 
                          Bin=unique(rates.bin[,"binAssignment"]), 
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
    bin.num  <- as.numeric(RECOM.BIN.MX[,"Bin"])
    ind.and  <- which(persist.mx[,i]%in%bin.num & persist.mx[,j]%in%bin.num)
    ntis.and <- persist.ntis[ind.and] 
    
    #subset of PERSIST.MX$hits with recombination rates for i AND j bins
    persist.and.mx <- persist.mx[ind.and,]
    i.and          <- persist.and.mx[,i]
    j.and          <- persist.and.mx[,j]
    
    #assign the subset of i and j bins to recombination rates from RECOM.BIN.MX
    #query   <- i or j bins with recombination rates
    #subject <- bin from RECOM.BIN.MX (stored as bin.num variable), 
    #which has recombination data (etc.) for a particular bin
    olap.i.and <- WhichOverlap(start.query=i.and, 
                               end.query=i.and, 
                               space.query=rep("a", length(i.and)),
                               start.subject=bin.num, 
                               end.subject=bin.num, 
                               space.subject=rep("a", length(bin.num)),
                               maxgap=-1L, minoverlap=1L)
    olap.j.and <- WhichOverlap(start.query=j.and, 
                               end.query=j.and, 
                               space.query=rep("a", length(j.and)),
                               start.subject=bin.num, 
                               end.subject=bin.num, 
                               space.subject=rep("a", length(bin.num)),
                               maxgap=-1L, minoverlap=1L)
    
    hits.i <- RECOM.BIN.MX[olap.i.and[,"subject"],]
    hits.j <- RECOM.BIN.MX[olap.j.and[,"subject"],]
    ij.and.mx <- cbind(Chr=hits.i[,"Chr"], 
                       i=i.and[olap.i.and[,"query"]], hits.i[,-c(1,2)], 
                       j=j.and[olap.j.and[,"query"]], hits.j[,-c(1,2)])
    
    #change column names
    ind.count                          <- which(colnames(ij.and.mx)=="COUNTINDVrb")
    colnames(ij.and.mx)[ind.count]     <- c("iCOUNTINDVrb", "jCOUNTINDVrb")
    
    ind.mean                           <- which(colnames(ij.and.mx)=="MEANrb")
    colnames(ij.and.mx)[ind.mean]      <- c("iMEANrb", "jMEANrb")
    
    ind.med                            <- which(colnames(ij.and.mx)=="MEDIANrb")
    colnames(ij.and.mx)[ind.med]       <- c("iMEDIANrb", "jMEDIANrb")
    
    ind.sdev                           <- which(colnames(ij.and.mx)=="SDEVrb")
    colnames(ij.and.mx)[ind.sdev]      <- c("iSDEVrb", "jSDEVrb")
    
    ind.indvrates                      <- which(colnames(ij.and.mx)=="INDVRATEScMperMb")
    colnames(ij.and.mx)[ind.indvrates] <- c("iINDVRATEScMperMb", "jINDVRATEScMperMb")
    
    ind.indvmap                        <- which(colnames(ij.and.mx)=="INDVMAPcM")
    colnames(ij.and.mx)[ind.indvmap]   <- c("iINDVMAPcM", "jINDVMAPcM")
    
    #statistical analyses pooling individual rates of i and j bins
    mx1.list <- list()
    mx1.list <- list(INDVRATEScMperMb=ij.and.mx[,ind.indvrates])
    result.indv  <- lapply(mx1.list, function(x){
      ij.pool        <- apply(x, MARGIN=1, function(x) paste(x, collapse=","))
      ij.pool        <- lapply(ij.pool, function(x) unlist(strsplit(x, split=",")) )
      ij.pool        <- lapply(ij.pool, function(x) as.numeric(x))
      mean.ij.pool   <- lapply(ij.pool, function(x) mean(x))
      median.ij.pool <- lapply(ij.pool, function(x) median(x))
      sd.ij.pool     <- lapply(ij.pool, function(x) sd(x)) #NA when n=1
      min.ij.pool    <- lapply(ij.pool, function(x) min(x))
      max.ij.pool    <- lapply(ij.pool, function(x) max(x))
      #P in column names to indicate pool of individual values for i and j bins
      #i=ij.and.mx[,c("i")], iCOUNTINDVrb=ij.and.mx[, ind.count[1]],
      #j=ij.and.mx[,c("j")], jCOUNTINDVrb=ij.and.mx[, ind.count[2]],
      all.val <- cbind(MEANijP= mean.ij.pool , MEDIANijP=median.ij.pool, 
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
      all.val <- cbind(ij.and.mx, NTISij=ntis.and, MEANij=mean.ij, MEDIANij=median.ij, 
                       SDEVij=sd.ij, MINij=min.ij, MAXij=max.ij, DIFFij=diff.ij, 
                       result.indv$INDVRATEScMperMb)
    })
    
    allMEANrb[[paste0("chr", chr)]]   <- result$MEANrb 
    allMEDIANrb[[paste0("chr", chr)]] <- result$MEDIANrb
    
    #incl  <- count of included recombination rates assigned to bins contacting OR
    #         count of ij contacts with corresponding rates
    #total <- total number of recombination rates OR ij contacts
    #perc  <- percentage
    incl  <- c()
    total <- c()
    perc  <- c()
    rates.stat <- c(incl=nrow(olap.mx), tot=nrow(recomRates.mx), 
                    perc=(nrow(olap.mx)/nrow(recomRates.mx))*100)
    ij.stat    <- c(incl=nrow(ij.and.mx), tot=nrow(persist.mx), 
                    perc=(nrow(ij.and.mx)/nrow(persist.mx)*100))
    
    allrates.stat[[paste0("chr", chr)]] <- rates.stat
    allij.stat[[paste0("chr", chr)]]    <- ij.stat
    
    #combine important results in a list 
    RECOM.PERSIST <- list(RECOM.PERSIST.MX=result, RECOM.BIN.MX=RECOM.BIN.MX,
                          RATESstat=round(rates.stat, digits=2), 
                          CONTACTSstat=round(ij.stat, digits=2))
    
    save(RECOM.PERSIST, file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_RecomVsPersist.RData"))
    
  } #end bracket for chr foreach
  
  allresult <- lapply(list(MEANrb=allMEANrb, MEDIANrb=allMEDIANrb), 
    function(x){
      x <- do.call(rbind, x)
      x <- as.data.frame(x)
  })
  
  allstat   <- lapply(list(RATESstat=allrates.stat, CONTACTSstat=allij.stat), 
                      function(x) do.call(rbind, x))
 
  allstat1  <- lapply(allstat, function(x){
    incl    <- sum(x[,"incl"])
    tot     <- sum(x[,"tot"])
    perc    <- (incl/tot)*100
    all.val <- cbind(incl, tot, perc)
  })
  
  allrates.stat <- round( as.vector(allstat1[[1]]), digits=2 )
  names(allrates.stat) <- c("incl", "tot", "perc")
  
  allij.stat    <- round( as.vector(allstat1[[2]]), digits=2 )
  names(allij.stat) <- c("incl", "tot", "perc")
  
  allstat       <- list(RATESstat=rbind(allstat$RATESstat, chrALL=allrates.stat), 
                         CONTACTSstat=rbind(allstat$CONTACTSstat, chrALL=allij.stat))
  
  RECOM.PERSIST <- list(RECOM.PERSIST.MX=allresult, STATpChr=allstat,
                        RATESstat=allrates.stat, CONTACTSstat=allij.stat)
  
  save(RECOM.PERSIST, file=paste0(output.dir, "/chrALL_min", gc, "Mb_RecomVsPersist.RData"))

} #end bracket for gc foreach

#rm(list=ls())

