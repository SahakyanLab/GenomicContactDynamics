#BaseAgeVsPersist
start_time <- Sys.time()
################################################################################
#Directories

#functions to source
#lib <- "/Users/ltamon/ProjectBoard/lib"
lib <- "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for task
#objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
objective.dir <- "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#Base age data 
#baseAge.dir <- "/Users/ltamon/Database/BaseAgeEns76_splitData"
baseAge.dir <- "/t1-home/icbml/ltamon/Database/BaseAgeEns76_splitData"
#HiC_Human21 persistent contact files
#persist.dir <- "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
persist.dir <- "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
#dir.create(paste0(objective.dir, "/out_BaseAgeVsPersist"))
output.dir <- paste0(objective.dir, "/out_BaseAgeVsPersist") 

################################################################################
#Packages

##CRAN
library(foreach)
library(data.table) #for fread(), rbindlist()
##BIOCONDUCTOR
library(IRanges) #for WhichOverlap

################################################################################
#Set values 

bin.length <- 40000L
#2(2Mb) or "05"(0-5 Mb), minimum gap acceptedbetween contacting bins 
#they should be far enough to filter for contacts within a TAD
#gc.v <- c(2, "05") 
gc <- 2
#chromosomes
#chr.v <- c(1:22, "X") 
chr <- 1
#column names of contacting bins in the persistent contact matrices
i <- "i"
j <- "j"

################################################################################ 
#Functions

source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))

################################################################################
################################################################################
#foreach(gc=gc.v, .inorder=TRUE) %do% {

  #allIJ.AND.MX   <- list()
  #allPOOL.AND.MX <- list()
  #allPOOL.EI.MX  <- list()

  #foreach(chr=chr.v, .inorder=TRUE) %do% {
      
    #PARTI - Assign base age to bins (only bins involved in contacts)
    
    baseAge.mx <- fread(file=paste0(baseAge.dir, "/chr", chr,"_BaseAgeEns76.txt"), 
                        header=FALSE, data.table=FALSE, drop=c(1,2))
    colnames(baseAge.mx) <- c("base", "ancestor", "score", "colourRGB")
    #no missing values in all base age datasets chr 1:2, X, Y and MT
    
    load(paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))
    persist.mx   <- PERSIST.MX$hits[,c(i, j)] #has the ij contacts
    persist.ntis <- PERSIST.MX$ntis #corresponding persistence score
    
    bins.uniq <- unique(unlist(persist.mx, use.names=FALSE))
    bins.uniq <- sort(as.numeric(bins.uniq))
    bin.end   <- bins.uniq*bin.length
    
    olap.mx <- WhichOverlap(start.query=baseAge.mx[,"base"], 
                            end.query=baseAge.mx[,"base"], 
                            space.query=rep("a", nrow(baseAge.mx)),
                            start.subject=bin.end-bin.length+1, 
                            end.subject=bin.end, 
                            space.subject=rep("a", length(bin.end)),
                            maxgap=0L, minoverlap=1L) 
    
    bin.age       <- cbind(bin=bins.uniq[olap.mx[,"subject"]], 
                        baseAge.mx[olap.mx[,"query"],])
    indv.base.b   <- by(bin.age[,"base"], bin.age[,"bin"], 
                        function(x) paste(x, collapse = ";"))
    indv.ancstr.b <- by(bin.age[,"ancestor"], bin.age[,"bin"], 
                        function(x) paste(x, collapse = ";"))
    indv.score.b  <- by(bin.age[,"score"], bin.age[,"bin"], 
                        function(x) paste(x, collapse = ";"))
    indv.RGB.b    <- by(bin.age[,"colourRGB"], bin.age[,"bin"], 
                        function(x) paste(x, collapse = ";"))
    
    BIN.AGE.MX <- cbind(bin=unique(bin.age[,"bin"]), 
                        base=indv.base.b,
                        ancestor=indv.ancstr.b,
                        score=indv.score.b,
                        colourRGB=indv.RGB.b)
    #---------------------------------------
    #PARTII - Intersect ij contacts and base age matrix  
    #take contacts that have baseAge data for (A) both i and j or (B) either i and j
    
    bin.num <- as.numeric(BIN.AGE.MX[,"bin"])
    ind.and <- which(persist.mx[,i]%in%bin.num & persist.mx[,j]%in%bin.num)
    ind.or  <- which(persist.mx[,i]%in%bin.num | persist.mx[,j]%in%bin.num)
    ind.ei  <- setdiff(ind.or, ind.and)
    #nme <- paste0("chr", chr)
  
    if(length(ind.or)==0){
      stop("All contacts have no base age data.")
    } else {
      #-----------------------
      #(A)i AND j have base age data 
      if(length(ind.and)==0){
        IJ.AND.MX   <- NULL
        POOL.AND.MX <- NULL
      } else {
        ntis.and <- persist.ntis[ind.and] 
        i.and    <- persist.mx[ind.and,i]
        j.and    <- persist.mx[ind.and,j]
        
        olap.i.and <- WhichOverlap(start.query=i.and, 
                                   end.query=i.and, 
                                   space.query=rep("a", length(i.and)),
                                   start.subject=bin.num, 
                                   end.subject=bin.num, 
                                   space.subject=rep("a", length(bin.num)),
                                   maxgap=0L, minoverlap=1L)
        olap.j.and <- WhichOverlap(start.query=j.and, 
                                   end.query=j.and, 
                                   space.query=rep("a", length(j.and)),
                                   start.subject=bin.num, 
                                   end.subject=bin.num, 
                                   space.subject=rep("a", length(bin.num)),
                                   maxgap=0L, minoverlap=1L)
        if(
          (length(i.and)+length(ntis.and)+
           nrow(olap.i.and)+nrow(olap.i.and))!=(length(i.and)*4)
          ){
            stop("Error in checkpoint comparing lengths.")
          }
  
        IJ.AND.MX <- cbind(BIN.AGE.MX[olap.i.and[,"subject"],], 
                           BIN.AGE.MX[olap.j.and[,"subject"],])
        
        #selecting columns to be pooled
        ind.base   <- which(colnames(IJ.AND.MX)=="base")
        ind.ancstr <- which(colnames(IJ.AND.MX)=="ancestor")
        ind.score  <- which(colnames(IJ.AND.MX)=="score")
        ind.RGB    <- which(colnames(IJ.AND.MX)=="colourRGB")
        
        #pool individual base age data for each contact
        mx.list <- list()
        #list of 2 column dataframes
        mx.list <- list(base=IJ.AND.MX[,ind.base], 
                        ancestor=IJ.AND.MX[,ind.ancstr],
                        score=IJ.AND.MX[,ind.score], 
                        colourRGB=IJ.AND.MX[,ind.RGB])
        pool.and <- lapply(mx.list, function(x){
          #returns vector of pooled strings from i and j
          pool <- apply(x, MARGIN=1, function(x) paste(x, collapse=";"))
        })
        
        POOL.AND.MX <- do.call(cbind, list(i=i.and[olap.i.and[,"query"]], 
                                           j=j.and[olap.j.and[,"query"]],
                                           ntis=ntis.and,
                                           do.call(cbind, pool.and)))
      } #else of ind.and if-else statement
      
      #-----------------------
      #(B)either i or j has base age data
      
      if(length(ind.ei)==0){
        POOL.EI.MX <- NULL
      } else {
        ntis.ei <- persist.ntis[ind.ei] 
        contacts.ei <- persist.mx[ind.ei,c(i,j)]
        lst <- list()
        for (a in 1:length(ind.ei)){
          cntct <- contacts.ei[a,]
          bin.withAge <- cntct[which(cntct%in%bin.num)][[1]]
          lst[[a]] <- BIN.AGE.MX[which(BIN.AGE.MX[,"bin"]==bin.withAge),]
        }
      
        POOL.EI.MX <- cbind(contacts.ei, ntis=ntis.ei, do.call(rbind, lst))
      } #else of ind.ei if-else statement
      
      #allIJ.AND.MX[[nme]]   <- IJ.AND.MX
      #allPOOL.AND.MX[[nme]] <- POOL.AND.MX
      #allPOOL.EI.MX[[nme]]  <- POOL.EI.MX
      
      BASEAGE.PERSIST <- list(BIN.AGE.MX=BIN.AGE.MX, 
                              POOL.AND.MX=POOL.AND.MX, 
                              POOL.EITHER.MX=POOL.EI.MX)
      
      save(BASEAGE.PERSIST, file=paste0(output.dir, "/chr", chr, "_min", gc, 
                                        "Mb_BaseAgeVsPersist.RData"))
    
    } #else of ind.or if-else statement
  
#  } #end bracket for chr foreach 
  
  #all_mx.lst <- lapply(list(allIJ.AND.MX=allIJ.AND.MX, 
  #                       allPOOL.AND.MX=allPOOL.AND.MX, 
  #                       allPOOL.EI.MX=allPOOL.EI.MX), 
  #             function(x){
  #               all_mx <- lapply(as.list(names(x)), function(y){
  #                 mx <- x[[y]]
  #                 cbind( chr=rep(y, nrow(mx)), mx )
  #               })
  #               do.call(rbind, all_mx)
  #             })
  
  #BASEAGE.PERSIST <- list(IJ.AND.MX=all_mx.lst$allIJ.AND.MX, 
  #                        POOL.AND.MX=all_mx.lst$allPOOL.AND.MX, 
  #                        POOL.EITHER.MX=all_mx.lst$allPOOL.EI.MX)
  
  #save(BASEAGE.PERSIST, file=paste0(output.dir, "/chrALL_min", gc, 
  #                                  "Mb_BaseAgeVsPersist.RData"))
#} #end bracket for gc foreach

end_time <- Sys.time()
end_time-start_time    

#rm(list=ls())
