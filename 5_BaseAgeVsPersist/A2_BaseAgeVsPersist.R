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
dir.create(paste0(objective.dir, "/out_BaseAgeVsPersist_new"))
output.dir = paste0(objective.dir, "/out_BaseAgeVsPersist_new") 


################################################################################
#Packages

##CRAN
library(foreach)
library(data.table) #for fread(), rbindlist()
##BIOCONDUCTOR
library(IRanges) #for WhichOverlap

################################################################################
#Set values 

bin.length = 40000L

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gc.v = "2" #c("2", "05") 

#chromosomes
chr.v = c(19, 21, 22) #c(1:22, "X") 

#gc="2"
#chr=21

################################################################################ 
################################################################################
foreach(chr=chr.v, .inorder=TRUE) %do% {
  baseAge.mx <- fread(file=paste0(baseAge.dir, "/chr", chr,"_BaseAgeEns76.txt"), 
                      header=FALSE, data.table=FALSE, drop=c(1,2))
  colnames(baseAge.mx) <- c("base", "ancestor", "score", "colourRGB")
  #no missing values in all base age datasets chr 1:2, X, Y and MT
  #actual base age is the second position
  
  for(gc in gc.v){

    load(paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))
    
    #get the relevant bins (from persist matrix) for that chromosome
    #bins sorted to increasing values so the ranges are also increasing
    bins.uniq <- unique( c(unique(PERSIST.MX$hits[,"i"]),
                           unique(PERSIST.MX$hits[,"j"])) )
    bins.uniq <- sort(bins.uniq, decreasing=FALSE)
    binOfAge <- ceiling(baseAge.mx[,"base"]/bin.length)
    
    #identify ages that fall in a contacting bin
    ageOnContact <- match(binOfAge, bins.uniq)
    ind.ageOnContact <- which(!is.na(ageOnContact))
    mx.age.bin <- cbind(bin=binOfAge[ind.ageOnContact], 
                        baseAge.mx[ind.ageOnContact,])
  
    count.indv.base.b <- by(mx.age.bin[,"base"], mx.age.bin[,"bin"], 
                            function(x) length(x))
    indv.base.b   <- by(mx.age.bin[,"base"], mx.age.bin[,"bin"], 
                        function(x) paste(x, collapse = ";"))
    indv.ancstr.b <- by(mx.age.bin[,"ancestor"], mx.age.bin[,"bin"], 
                        function(x) paste(x, collapse = ";"))
    indv.score.b  <- by(mx.age.bin[,"score"], mx.age.bin[,"bin"], 
                        function(x) paste(x, collapse = ";"))
    indv.RGB.b    <- by(mx.age.bin[,"colourRGB"], mx.age.bin[,"bin"], 
                        function(x) paste(x, collapse = ";"))
    
    bin.assigned <- as.numeric(names(count.indv.base.b))
    
    if(identical(bin.assigned, unique(mx.age.bin[,"bin"]))==FALSE){
      stop("Error in checkpoint before making AGE.BIN.MX.")
    }
    
    AGE.BIN.MX <- cbind(bin=bin.assigned,
                        count=count.indv.base.b,
                        base=indv.base.b,
                        ancestor=indv.ancstr.b,
                        score=indv.score.b,
                        colourRGB=indv.RGB.b)
    
    #---------------------------------------
    #PARTII - Intersect ij contacts and base age matrix  
    #take contacts that have baseAge data for (A) both i and j or (B) either i and j
    
    ##*** eliminate persist.mx, use the original object
    
    ind.and <- which(PERSIST.MX$hits[,"i"]%in%bin.assigned & PERSIST.MX$hits[,"j"]%in%bin.assigned)
    ind.or <- which(PERSIST.MX$hits[,"i"]%in%bin.assigned | PERSIST.MX$hits[,"j"]%in%bin.assigned)
    ind.ei  <- setdiff(ind.or, ind.and)

    if(length(ind.or)==0){
      stop("All contacts have no base age data.")
    } else {
      #-----------------------
      #(A)i AND j have base age data 
      if(length(ind.and)==0){
        POOL.AND.MX <- NULL
      } else {

        #bin.assigned same order as AGE.BIN.MX
        i.match <- match(PERSIST.MX$hits[ind.and,"i"], bin.assigned)
        j.match <- match(PERSIST.MX$hits[ind.and,"j"], bin.assigned)
        
        mx.list <- list()
        columns <- c("count", "base", "ancestor", "score", "colourRGB")
        for(col in columns){
          mx <- cbind(AGE.BIN.MX[i.match, col], 
                      AGE.BIN.MX[j.match, col])
          mx.list[[col]] <- apply(mx, MARGIN=1, 
                                  function(x) paste(x, collapse=";"))
        }
        
        HiCsignalStr.and <- apply(PERSIST.MX$hits[ind.and,3:23], MARGIN=1, 
                                  function(x){paste(x, collapse=";")})
        POOL.AND.MX <- cbind( i=PERSIST.MX$hits[ind.and,"i"],
                              j=PERSIST.MX$hits[ind.and,"j"],
                              do.call(cbind, mx.list), 
                              ntis=PERSIST.MX$ntis[ind.and],
                              valsum=PERSIST.MX$valsum[ind.and],
                              Co_Hi_Lu_LV_RV_Ao_PM_Pa_Sp_Li_SB_AG_Ov_Bl_MesC_MSC_NPC_TLC_ESC_FC_LC=HiCsignalStr.and )
      
      } #else of ind.and if-else statement
      
      #-----------------------
      #(B)either i or j has base age data
      
      if(length(ind.ei)==0){
        POOL.EI.MX <- NULL
      } else {
        contacts.ei <- PERSIST.MX$hits[ind.ei,c("i", "j")]
        lst <- list()
        for (a in 1:length(ind.ei)){
          cntct <- contacts.ei[a,]
          bin.withAge <- cntct[which(cntct%in%bin.assigned)][[1]]
          lst[[a]] <- AGE.BIN.MX[which(AGE.BIN.MX[,"bin"]==bin.withAge),]
        }
        
        HiCsignalStr.ei <- apply(PERSIST.MX$hits[ind.ei,3:23], MARGIN=1, 
                                  function(x){paste(x, collapse=";")})
        POOL.EI.MX <- cbind(contacts.ei,
                            do.call(rbind, lst),
                            ntis=PERSIST.MX$ntis[ind.ei],
                            valsum=PERSIST.MX$valsum[ind.ei],
                            Co_Hi_Lu_LV_RV_Ao_PM_Pa_Sp_Li_SB_AG_Ov_Bl_MesC_MSC_NPC_TLC_ESC_FC_LC=HiCsignalStr.ei )
        
      } #else of ind.ei if-else statement
      
      BASEAGE.PERSIST <- list(AGE.BIN.MX=AGE.BIN.MX, 
                              POOL.AND.MX=POOL.AND.MX, 
                              POOL.EITHER.MX=POOL.EI.MX)
      save(BASEAGE.PERSIST, file=paste0(output.dir, "/min", gc, "Mb_chr", chr,  
                                        "_BaseAgeVsPersist.RData"))
    } #else of ind.or if-else statement
  }
}

end_time <- Sys.time()
end_time-start_time    

#rm(list=ls())