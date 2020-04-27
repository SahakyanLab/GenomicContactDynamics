################################################################################
# Map any feature to HiC contact persist bins
# Parallel execution by chromosome
# DEPENDENCIES:
## library(data.table)
## library(foreach)
## library(itertools)
## library(doParallel)
## source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
## source(paste0(lib, "/loadRData.R"))
## source(paste0(lib, "/LO_mapToHiCcontactPersistBins.R")) 
### FUNCTION ###################################################################
mapToHiCcontactPersistBins <- function(
  
  # PERSIST.MX object directory
  persist.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc",
  
  ## Feature directory
  feat.dir = feat.dir,
  # Separate files per chromosome?
  featSepByChr = FALSE,
  #---- if featSepByChr = TRUE
  # Program will look for each file on feat.dir using chr number eg. *chrX* as pattern
  # Beware of duplicate files for each chr
  
  #---- if featSepByChr = FALSE
  # Path to file
  featurefile = paste0(feat.dir, "/hg19_TOP2B_MCF7_GSE66753_atleast_two_peaks.bed"),
  # Combine output of chromosomes?
  combineOut = TRUE,
  
  feat.header = TRUE,
  # Column numbers or names
  start.coord = "start",
  end.coord = "end",
  # Table coordinate system
  # zero-based (0-start, half-open) or one-based (1-start, fully-closed (1-based))
  # Coordinates stored in UCSC Genome Browser tables are zero-based while those on 
  # web interfaces are one-based
  feat.coordSys = "one-based",
  chr.col = "chr",
  # Remove irrelevant columns
  # Recommended that the featurefile should have a column uniquely identifying each entry/row
  drop.col = NULL, 
  # Column names after dropping chr.col and drop.col (follow order in original table)
  # If NULL, keep original column names
  coluNames = NULL,
  
  # Output directory
  output.dir = paste0(objective.dir, "/out_mTPB_topo"),
  # Output identifier
  out.name = "MCF7_TOP2B_within_40KbBin",
  
  # HiC resolution
  bin.len = 40000L,
  # 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap between contacting bins
  gcb = "min2Mb", 
  chr.v = paste("chr", c(1:22, "X"), sep=""),
  nCPU = 23L,
  
  # Overlap parameters
  min.olap = 1L,
  max.gap = -1L,
  type.olap = "any",
  
  # LiftOver
  doLiftOver = TRUE,
  LOchain = "hg19ToHg38"
){

  # Initialise the nCPU-driven (number of CPUs) setting for do/dopar (foreach)
  if(nCPU > 1){
    registerDoParallel(cores=nCPU)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  chr.v.len <- length(chr.v)
  
  if( featSepByChr==FALSE & grepl(pattern=".RData", x=featurefile) ){
    feat.df <- loadRData(fileName=featurefile)
  } else {
    feat.df <- fread(file=featurefile, 
                     header=feat.header, data.table=FALSE, 
                     stringsAsFactors=FALSE)
  }
  
  toExport <- c(
    "objective.dir",
    "persist.dir",
    "bin.len",
    
    "feat.dir",
    "featSepByChr",
    "combineOut",
    "feat.header",
    "feat.df",
    "start.coord",
    "end.coord",
    "chr.col",
    "drop.col",
    "coluNames",
    
    "output.dir",
    "out.name",
    "max.gap",
    "min.olap",
    "type.olap",
    "chr.v",
    "gcb"
  )
  
  #### PARALLEL EXECUTION #########
  
  FEATURE.BIN <- foreach( itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                          .inorder=TRUE, .combine="rbind",
                          .export=toExport, 
                          .noexport=ls()[!ls()%in%toExport]
                             
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
    
      chr <- chr.v[i]
      
      if(featSepByChr==TRUE){
        
        featurefile=paste0(feat.dir, "/", dir(feat.dir, pattern=chr))
        
        if( grepl(".RData", x=featurefile) ){
          feat.df <- loadRData(fileName=featurefile)
        } else {
          feat.df <- fread(file=featurefile, header=feat.header, 
                           data.table=FALSE, stringsAsFactors=FALSE)
        }
        
      } else {
        feat.df <- feat.df[feat.df[[chr.col]]==chr,]
      }
      
      # Check for missing values in feat.df
      colNASums <- colSums(x=is.na(feat.df))
      if(any(colNASums!=0)){
        print( paste0(chr, ":Missing value/s in column/s ", 
                      paste(which(colNASums!=0), collapse=","), ".") )
      }
      
      # Change column names of feature coordinates/position
      if(start.coord==end.coord){
        
        #setnames(feat.df, old=colnames(feat.df[,c(start.coord, chr.col)])[1], 
        #         new="pos")
        
        # Start coord can be column name or number
        setnames(x=feat.df, old=start.coord, new="pos")
        
      } else {
        
        # Lengths of features
        size <- feat.df[[end.coord]]-feat.df[[start.coord]]
        # Just checking if end always greater than start as assummed 
        if( sum(size<0) ){
          stop(paste0("Chr", chr,
                      ":There are start coordinates > end coordinates."))
        } 
        
        #setnames(feat.df, old=colnames(feat.df[,c(start.coord, chr.col)])[1], 
        #         new="start")
        setnames(x=feat.df, old=start.coord, new="start")
        #setnames(feat.df, old=colnames(feat.df[,c(end.coord, chr.col)])[1], 
        #         new="end")
        setnames(x=feat.df, old=end.coord, new="end")
        
      }
      
      # Remove irrelevant columns
      if( is.numeric( c(chr.col, drop.col) ) ){
        feat.df <- feat.df[,-drop.ind]
      } else {
        drop.ind <- which(colnames(feat.df)%in%c(chr.col, drop.col))
        feat.df <- feat.df[,-drop.ind]
      }
      
      load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      # Get the relevant bins for that chromosome (from persist matrix) 
      bins.uniq <- unique( c(unique(PERSIST.MX$hits[,"i"]),
                             unique(PERSIST.MX$hits[,"j"])) )
      
      rm("PERSIST.MX"); gc()
      
      # Assign features to bins of the chromosome
      bin.end   <- bins.uniq*bin.len
      bin.start <- bin.end-bin.len+1
      
      NAbin <- NULL
      
      if(doLiftOver){
        # Corresponds to order of bin.end, bin.start hence bins.uniq
        LO.mx <- LO_mapToHiCcontactPersistBins(out.dir=out.dir, 
                                               out.name=paste0(chr, "_Persist_", gcb),
                                               chr=chr,
                                               start=bin.start, end=bin.end,
                                               LOchain=LOchain, rmchain=FALSE)
        NAbin.TF <- is.na(LO.mx[,"end"]) & is.na(LO.mx[,"start"]) 
        bin.end <- LO.mx[!NAbin.TF,"end"] 
        bin.start <- LO.mx[!NAbin.TF,"start" ] 
        NAbin <- bins.uniq[NAbin.TF]
        bins.uniq <- bins.uniq[!NAbin.TF]; rm(NAbin.TF)
        
        if( unique(c(length(bin.end), length(bin.start))!=length(bins.uniq)) ){
          stop("Checkpoint: LiftOver output missing bins.")
        }
        rm(LO.mx); gc()
      }
      
      if(start.coord==end.coord){
        
        feat.start <- feat.end <- feat.df$pos
        
      } else {
        
        # Change start coordinates according to feature coordinate system 
        if(feat.coordSys=="zero-based"){
          feat.start <- (feat.df$start) + 1L
          print("Feature file: 0-based coordinate system, start coordinate + 1.")
        } else if(feat.coordSys=="one-based"){
          feat.start <- feat.df$start
          print("Feature file: 1-based coordinate system, start coordinate maintained.")
        } else {
          stop("Invalid input for feature coordinate system.")
        }
        
        feat.end <- feat.df$end
        
      }
      
      # Query   <- features (each row of feat.df)
      # Subject <- persist bins (bins.uniq)
      olap.mx <- WhichOverlap(start.query=feat.start, 
                              end.query=feat.end, 
                              space.query=rep("a", length( feat.df[[1]] )),
                              start.subject=bin.start, 
                              end.subject=bin.end, 
                              space.subject=rep("a",length(bin.end)),
                              maxgap=max.gap, minoverlap=min.olap,
                              type=type.olap)
      
      colString.lst <- sapply(X=colnames(feat.df), simplify=FALSE, 
                              FUN=function(col){
                                by(data=feat.df[olap.mx[,"query"], col],
                                   INDICES=bins.uniq[ olap.mx[,"subject"] ],
                                   FUN=function(x) paste(x, collapse = ";"))
                              })
      
      df <- do.call("cbind", colString.lst)
      rm("colString.lst", "bin.end", "bin.start"); gc()
      
      rownames(df) <- NULL
      if( !is.null(coluNames) ){
        colnames(df) <- coluNames
      }
      
      # Count data per bin
      countPerBin <- by(data=olap.mx[,"query"],
                        INDICES=bins.uniq[ olap.mx[,"subject"] ],
                        # Unique is optional because query-subj pair in olap.mx
                        # should be unique always
                        FUN=function(x) length(unique(x)))
      
      rm("olap.mx"); gc()
      
      FEATURE.BIN.MX <- cbind.data.frame(chr=rep(chr),
                                         bin=as.numeric(names(countPerBin)),
                                         df,
                                         countPerBin=as.numeric(countPerBin),
                                         stringsAsFactors=FALSE)
      if( length(NAbin)!=0 ){
        appendd <- data.frame(matrix(data=NA, nrow=length(NAbin), ncol=ncol(FEATURE.BIN.MX),
                                     dimnames=list(NULL, c("chr", "bin", colnames(df), "countPerBin"))
                                     )
                              )
        appendd$bin <- NAbin
        appendd$chr <- chr
        
        FEATURE.BIN.MX <- rbind(FEATURE.BIN.MX, appendd)
        rm(appendd, NAbin)
      }
      
      rm(df, countPerBin); gc()
      
      if(combineOut==FALSE){
        save(FEATURE.BIN.MX, file=paste0(output.dir, "/", chr, "_", 
                                         gcb, "_", out.name, ".RData"))
      } else {
        return(FEATURE.BIN.MX)
      }
      
      print(paste0(chr, " done!"))
      
    }) # itr sapply end
    
    return(do.call(rbind, chunk))
    
  } # foreach end
  
  ### END OF PARALLEL EXECUTION ###
  
  if(combineOut){
    FEATURE.BIN.MX <- FEATURE.BIN; rm(FEATURE.BIN); gc()
    save(FEATURE.BIN.MX, file=paste0(output.dir, "/chrALL_", 
                                     gcb, "_", out.name, ".RData"))
  } else {
    rm(FEATURE.BIN); gc()
  }
  
} # Function end
################################################################################
suppressWarnings(suppressPackageStartupMessages(library(compiler)))
mapToHiCcontactPersistBins <- cmpfun(mapToHiCcontactPersistBins, options=list(suppressUndefined=TRUE))
################################################################################