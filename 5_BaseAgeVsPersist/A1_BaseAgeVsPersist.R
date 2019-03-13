# BaseAgeVsPersist
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
                            # "AlexMac", "AlexCluster"

#*** See whether the content of the <lib>, specified below, can be moved to the
#*** actual project folder for self consistency. Then, specify the location of
#*** <lib> within the whorunsit block below. If you put the dependent function
#*** in <lib> into the <objective.dir>/lib location, then the following line
#*** outside the whorunsit block, after output.dir specification is enough:
#***    lib = paste0(objective.dir, "/lib")

# Functions to source
# lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelCluster"){
    # Main directory for the task
    objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
    # Base age data
    baseAge.dir   = "/t1-home/icbml/ltamon/Database/BaseAgeEns76_splitData"
    # HiC_Human21 persistent contact files
    persist.dir   = "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelMac"){
    objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
    baseAge.dir   = "/Users/ltamon/Database/BaseAgeEns76_splitData"
    persist.dir   = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "AlexMac"){
    objective.dir = "./"
    baseAge.dir   = "/Volumes/Data/Database/BaseAgeEns76_splitData"
    persist.dir   = "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}

output.dir = paste0(objective.dir, "/out_BaseAgeVsPersist")
### OTHER SETTINGS #############################################################
# "2"(2MB gap) or "05"(0.5 MB minimum gap), refers to a minimum gap accepted to
# further process a given contact: two points should be far enough to filter for
# contacts within a TAD.
gc.v = c("2", "05")

# chromosomes
chr.v = c(1:22, "X")
bin.length = 40000L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(data.table) #for fread(), rbindlist()
library(IRanges)    #required for WhichOverlap()
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
start_time <- Sys.time()

foreach(chr=chr.v, .inorder=TRUE) %do% {
  baseAge.mx <- fread(file=paste0(baseAge.dir, "/chr", chr,"_BaseAgeEns76.txt"),
                      header=FALSE, data.table=FALSE, drop=c(1,2))
  colnames(baseAge.mx) <- c("base", "ancestor", "score", "colourRGB")
  binOfAge <- ceiling(baseAge.mx[,"base"]/bin.length)
  #no missing values in all base age datasets chr 1:2, X, Y and MT
  #actual base age is the second position

  #*** beware of gc() being garbage collector
  for(gc in gc.v){
    load(paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))

    #*** ~10.3 GB for chr1

    #*** For clarity, always add one space after <#> in your comments. This is
    #*** not required while commenting out lines of code, so only for comments.

    # get the relevant bins (from persist matrix) for that chromosome
    bins.uniq <- unique( c(unique(PERSIST.MX$hits[,"i"]),
                           unique(PERSIST.MX$hits[,"j"])) )

    #*** sorting does not mattter in the new way of the script, so depreciating the below line:
    # bins sorted to increasing values so the ranges are also increasing
    #bins.uniq <- sort(bins.uniq, decreasing=FALSE)
    #*** taking binOfAge operation out of the for loop since it does not depend in PERSIST.MX

    # identify ages that fall in a contacting bin
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
