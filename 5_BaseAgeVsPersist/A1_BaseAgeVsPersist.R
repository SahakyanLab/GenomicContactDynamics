#BaseAgeVsPersist
start_time <- Sys.time()

################################################################################
#Directories

#functions to source
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for task
#objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#Base age data 
#baseAge.dir = "/Users/ltamon/Database/BaseAgeEns76_splitData"
baseAge.dir = "/t1-home/icbml/ltamon/Database/BaseAgeEns76_splitData"
#HiC_Human21 persistent contact files
#persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
persist.dir = "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
#***** repetition below
#dir.create(paste0(objective.dir, "/out_BaseAgeVsPersist_sample"))
#output.dir = paste0(objective.dir, "/out_BaseAgeVsPersist_sample") 
output.dir = paste0(objective.dir, "/out_BaseAgeVsPersist") 

library(foreach)
library(data.table) #for fread(), rbindlist()
library(IRanges) #for WhichOverlap

bin.length = 40000L

#2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gc.v = c("2", "05") 

#chromosomes
chr.v = c(1:22, "X") 

#gc="2"
#chr=21

################################################################################ 
################################################################################
foreach(chr=chr.v, .inorder=TRUE) %do% {
  baseAge.mx <- fread(file=paste0(baseAge.dir, "/chr", chr,"_BaseAgeEns76.txt"), 
                      header=FALSE, data.table=FALSE, drop=c(1,2))
  colnames(baseAge.mx) <- c("base", "ancestor", "score", "colourRGB")
  binOfAge <- ceiling(baseAge.mx[,"base"]/bin.length)
  #no missing values in all base age datasets chr 1:2, X, Y and MT
  #actual base age is the second position
  
  #**** beware of gc() being garbage collector
  for(gc in gc.v){
    load(paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))
    
    #*** ~10.3 GB for chr1
    
    #get the relevant bins (from persist matrix) for that chromosome
    #bins sorted to increasing values so the ranges are also increasing
    bins.uniq <- unique( c(unique(PERSIST.MX$hits[,"i"]),
                           unique(PERSIST.MX$hits[,"j"])) )
    
    #*** sorting does not mattter in the new way of the script, so depreciating the below line:
    #bins.uniq <- sort(bins.uniq, decreasing=FALSE)
    #**** taking binOfAge operation out of the for loop since it does not depend in PERSIST.MX
    
    #identify ages that fall in a contacting bin
    #*** skipping the  creation of a PERMANENT ageOnContact object:
    #ageOnContact <- match(binOfAge, bins.uniq)
    ind.ageOnContact <- which(!is.na(  match(binOfAge, bins.uniq)  ))
    
    #*** No need for a redundant matrix, as the components are already there
    #mx.age.bin <- cbind(  bin=binOfAge[ind.ageOnContact], 
    #                      baseAge.mx[ind.ageOnContact,]   )
  
    #*** restructuring the code to skip the usage of mx.age.bin  
    count.indv.base.b <- by( data = baseAge.mx[ind.ageOnContact,"base"],
                             INDICES = binOfAge[ind.ageOnContact], 
                             FUN = function(x) length(x) )

    indv.base.b   <- by( data = baseAge.mx[ind.ageOnContact,"base"],
                         INDICES = binOfAge[ind.ageOnContact], 
                         FUN = function(x) paste(x, collapse = ";") )
    
    indv.ancstr.b <- by( data = baseAge.mx[ind.ageOnContact,"ancestor"],
                         INDICES = binOfAge[ind.ageOnContact], 
                         FUN = function(x) paste(x, collapse = ";") )
    
    indv.score.b  <- by( data = baseAge.mx[ind.ageOnContact,"score"],
                         INDICES = binOfAge[ind.ageOnContact], 
                         FUN = function(x) paste(x, collapse = ";") )
    
    indv.RGB.b    <- by( data = baseAge.mx[ind.ageOnContact,"colourRGB"],
                         INDICES = binOfAge[ind.ageOnContact], 
                         FUN = function(x) paste(x, collapse = ";") )
    bin.assigned <- as.numeric(names(count.indv.base.b))
    if( identical(bin.assigned, unique(binOfAge[ind.ageOnContact]) )==FALSE){
      stop("Error in checkpoint before making BASEAGE.BIN.MX.")
    }
    BASEAGE.BIN.MX <- cbind(bin=bin.assigned,
                            count=count.indv.base.b,
                            base=indv.base.b,
                            ancestor=indv.ancstr.b,
                            score=indv.score.b,
                            colourRGB=indv.RGB.b)
    #*** thats is a good idea to wrap up a part of the code with only a certain 
    # object creation, then deleting all the rest, with a non-compulsory
    # garbage collection (gc()):
    rm( list=c("count.indv.base.b",
               "indv.base.b", "indv.ancstr.b",
               "indv.score.b", "indv.RGB.b",
               "ind.ageOnContact", "bins.uniq") )
    save(BASEAGE.BIN.MX, file=paste0(output.dir, "/min", gc, "Mb_chr", chr,  
                                     "_bin40Kb_BaseAgeVsPersist.RData"))
  }
}

end_time <- Sys.time()
end_time-start_time    

#rm(list=ls())