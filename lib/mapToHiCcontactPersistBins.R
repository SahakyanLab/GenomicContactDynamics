################################################################################
# Map any feature to HiC contact persist bins; parallel execution by chromosome
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(data.table)
# library(foreach)
# library(itertools)
# library(doParallel)
# library(GenomicRanges)
# library(compiler)
# source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
# source(paste0(lib, "/loadRData.R"))
# source(paste0(lib, "/LO_mapToHiCcontactPersistBins.R")) 
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
  pattern.pref = NULL,
  pattern.suff = NULL,
  
  #---- if featSepByChr = FALSE
  # Path to file
  featurefile = paste0(feat.dir, "/hg19_TOP2B_MCF7_GSE66753_atleast_two_peaks.bed"),
  # Combine output of chromosomes?
  combineOut = FALSE,
  
  feat.header = NULL,
  # Column numbers or names
  start.coord = "start",
  end.coord = "end",
  # Table coordinate system
  # zero-based (0-start, half-open) or one-based (1-start, fully-closed (1-based))
  # Coordinates stored in UCSC Genome Browser tables are zero-based while those on 
  # web interfaces are one-based
  feat.coordSys = "zero-based",
  chr.col = "chr",
  # Remove irrelevant columns
  # Recommended that the featurefile should have a column uniquely identifying each entry/row
  drop.col = NULL, 
  # Column names after dropping chr.col and drop.col (follow order in original table)
  # If NULL, keep original column names
  coluNames = NULL,
  
  # Output directory
  out.dir = paste0(objective.dir, "/out_mTPB_topo"),
  # Output identifier
  out.name = "MCF7_TOP2B_within_40KbBin",
  
  # HiC resolution
  bin.len = 40000L,
  # 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap between contacting bins
  gcb = "min2Mb", 
  chr.v = paste("chr", c(1:22, "X"), sep=""),
  nCPU = 1L,
  
  # Overlap parameters
  min.olap = 1L,
  max.gap = -1L,
  type.olap = "any",
  
  # LiftOver
  doLiftOver = FALSE,
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
  
  if( !featSepByChr ){
    if(grepl(pattern=".RData", x=featurefile)){
      FEAT.df <- loadRData(fileName=featurefile)
    } else {
      FEAT.df <- fread(file=featurefile, 
                       header=feat.header, data.table=FALSE, 
                       stringsAsFactors=FALSE)
    }
    chr.v <- intersect( chr.v, unique(FEAT.df[[chr.col]]) )
    FEAT.df <- FEAT.df[FEAT.df[[chr.col]]%in%chr.v,]
  }
  
  toExport <- c(
    "objective.dir",
    "persist.dir",
    "bin.len",
    
    "feat.dir",
    "featSepByChr",
    "combineOut",
    "feat.header",
    "FEAT.df",
    "start.coord",
    "end.coord",
    "chr.col",
    "drop.col",
    "coluNames",
    
    "out.dir",
    "out.name",
    "max.gap",
    "min.olap",
    "type.olap",
    "chr.v",
    "gcb"
  )
  chr.v.len <- length(chr.v)
  
  #### PARALLEL EXECUTION #########
  FEATURE.BIN <- foreach( itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                          .inorder=TRUE, .combine="rbind",
                          .export=toExport, 
                          .noexport=ls()[!ls()%in%toExport]
                             
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
    
      chr <- chr.v[i]
      
      if(featSepByChr){
        featurefile=paste0(feat.dir, "/", dir(feat.dir, pattern=paste0(pattern.pref, chr, pattern.suff)))
        if( grepl(".RData", x=featurefile) ){
          feat.df <- loadRData(fileName=featurefile)
        } else {
          feat.df <- fread(file=featurefile, header=feat.header, 
                           data.table=FALSE, stringsAsFactors=FALSE)
        }
      } else {
        feat.df <- FEAT.df[FEAT.df[[chr.col]]==chr,]
      }
     
      if(unique(feat.df[[chr.col]])!=chr){ stop(paste0(chr, ":Feature ranges from different chr.")) }
      
      # Convert column numbers/indices to names, if applicable
      lapply(X=list(start.coord, end.coord, chr.col, drop.col), FUN=is.numeric)
      
      for(objct in c("start.coord", "end.coord", "chr.col", "drop.col") ){
        
        eval(parse(text=paste0( "is_num <- is.numeric(", objct, ")")))
        if(is_num){
          eval(parse(text=paste0( objct, " <- colnames(feat.df)[", objct, "]")))
        }
        
      }
      
      # Remove irrelevant columns
      feat.df <- feat.df[,setdiff(colnames(feat.df), c(chr.col, drop.col))]
      
      # Check for missing values in feat.df
      colNASums <- colSums(x=is.na(feat.df))
      if( any(colNASums>0) ){
        print( paste0(chr, ":Missing value/s in column/s ", 
                      paste(which(colNASums!=0), collapse=","), ".") )
      }
      
      # Change column names of feature coordinates/position
      if(start.coord==end.coord){
        setnames(x=feat.df, old=start.coord, new="pos")
      } else {
        # Lengths of features
        size <- feat.df[[end.coord]]-feat.df[[start.coord]]
        # Just checking if end always greater than start as assummed 
        if( any(size<0) ){
          stop(paste0("Chr", chr,
                      ":There are start coordinates > end coordinates."))
        }; rm(size)
        setnames(x=feat.df, old=start.coord, new="start")
        setnames(x=feat.df, old=end.coord, new="end")
      }
      
      load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      # Get the relevant bins for that chromosome (from persist matrix) 
      bins.uniq <- unique( c(unique(PERSIST.MX$hits[,"i"]),
                             unique(PERSIST.MX$hits[,"j"])) )
      
      rm(PERSIST.MX); gc()
      
      # Assign features to bins of the chromosome
      bin.end   <- bins.uniq*bin.len
      bin.start <- bin.end-bin.len+1
      
      NAbin <- NULL
      if( doLiftOver & !is.null(LOchain) ){
        print(paste0(chr, ":", LOchain, " liftover..."))
        # LO.mx$group corresponds to index of bin in bins.uniq
        LO.mx <- LO_mapToHiCcontactPersistBins(out.dir=out.dir, 
                                               out.name=paste0(chr, "_Persist_", gcb),
                                               chr=chr, bin.len=bin.len, bin.num=bins.uniq,
                                               start=bin.start, end=bin.end,
                                               LOchain=LOchain, rmchain=FALSE)
        # Bins not converted
        NAbin <- setdiff(bins.uniq, LO.mx$group)
        # bins.uniq can be duplicated when liftover of bin result into multiple ranges
        bins.uniq <- LO.mx$group 
        bin.start <- LO.mx$start
        bin.end <- LO.mx$end
        rm(LO.mx); gc()
      } else if( doLiftOver & is.null(LOchain) ){
        stop(paste0(chr, ":Provide LOchain."))
      } 
      
      if(start.coord==end.coord){
        feat.start <- feat.end <- feat.df$pos
      } else {
        # Change start coordinates according to feature coordinate system 
        if(feat.coordSys=="zero-based"){
          feat.start <- (feat.df$start) + 1L
          print(paste0(chr, ":Feature file: 0-based coordinate system, start coordinate + 1."))
        } else if(feat.coordSys=="one-based"){
          feat.start <- feat.df$start
          print(paste0(chr, ":Feature file: 1-based coordinate system, start coordinate maintained."))
        } else {
          stop(paste0(chr, ":Invalid input for feature coordinate system."))
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
                                   INDICES=bins.uniq[olap.mx[,"subject"]],
                                   FUN=function(x) paste(x, collapse = ";"))
                              })
      df <- do.call("cbind", colString.lst)
      rm(colString.lst, bin.end, bin.start, feat.df); gc()
      
      rownames(df) <- NULL
      if( !is.null(coluNames) ){
        colnames(df) <- coluNames
      }
      
      # Count data per bin
      countPerBin <- by(data=olap.mx[,"query"],
                        INDICES=bins.uniq[olap.mx[,"subject"]],
                        # Unique is optional because query-subj pair in olap.mx
                        # should be unique always
                        FUN=function(x) length(unique(x)))
      
      rm(olap.mx); gc()
      
      FEATURE.BIN.MX <- cbind.data.frame(chr=rep(chr),
                                         bin=as.numeric(names(countPerBin)),
                                         df,
                                         countPerBin=as.numeric(countPerBin),
                                         stringsAsFactors=FALSE)
      # NA bins from liftover plus bins with no overlaps
      NAbin <- unique( c(NAbin, setdiff(bins.uniq, FEATURE.BIN.MX$bin)) )
      if( length(NAbin)>0 ){
        appendd <- data.frame(matrix(data=NA, nrow=length(NAbin), ncol=ncol(FEATURE.BIN.MX),
                                     dimnames=list(NULL, c("chr", "bin", colnames(df), "countPerBin"))
                                     ), check.names=F 
                              )
        appendd$bin <- NAbin
        appendd$chr <- chr
        appendd$countPerBin <- 0
        FEATURE.BIN.MX <- rbind(FEATURE.BIN.MX, appendd)
        rm(appendd, NAbin)
      }
      if( any(duplicated(FEATURE.BIN.MX$bin)) ){ stop("Duplicated bins.") }
      rm(df, countPerBin); gc()
      
      FEATURE.BIN.MX <- FEATURE.BIN.MX[order(FEATURE.BIN.MX$bin, decreasing=FALSE),]
      
      print(paste0(chr, " done!"))
      
      if( combineOut==FALSE | is.null(combineOut) ){
        save(FEATURE.BIN.MX, file=paste0(out.dir, "/", chr, "_", 
                                         gcb, "_", out.name, ".RData"))
        return(NULL)
      } else if(combineOut){
        return(FEATURE.BIN.MX)
      } else {
        stop("combineOut: Invalid input")
      }
      
    }) # itr sapply end
    
    return(do.call(rbind, chunk))
    
  } # foreach end
  ### END OF PARALLEL EXECUTION ###
  
  if(combineOut){
    FEATURE.BIN.MX <- FEATURE.BIN; rm(FEATURE.BIN); gc()
    save(FEATURE.BIN.MX, file=paste0(out.dir, "/chrALL_", 
                                     gcb, "_", out.name, ".RData"))
  } else {
    rm(FEATURE.BIN); gc()
  }
  
} # Function end
################################################################################
mapToHiCcontactPersistBins <- cmpfun(mapToHiCcontactPersistBins, options=list(suppressUndefined=TRUE))
################################################################################