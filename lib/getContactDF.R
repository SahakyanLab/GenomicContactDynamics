################################################################################
# Get dataframe of upper triangle contacts (i-j column) including a column for 
# values based on metric and a column specifying whether the contact will be 
# included or not based on specified bins to be included/masked. Note that
# Cs.raw, Cs.norm, Cp, and SIM contacts can contain 0s. Cs=0 means that the
# contact is not in given tissue, Cp=0 means contact not in any tissue or not
# in given tissue. SIM=0 means contact not observed given cut-off for a contact.
# CII categorized = 0 means middle group. Function will return 0 values and
# will specify whether they are included or not based on contact filtering 
# arguments, but it's up to the downstream code how to deal with these.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(compiler)
# library(data.table)
# library(reshape2)
# source(paste0(lib, "/simulation_lib/filterContacts.R"))
# source(paste0(lib, "/categoriseValues.R"))
### FUNCTION ###################################################################
getContactDF <- function(metric.dir='', metric='', gcb = '', chr = '', 
                         ct = 'cell type', gap.range='closed vector',
                         incl.bin.x = NULL, incl.bin.y = NULL, # list, 
                         mask.bin.x = NULL, mask.bin.y = NULL,  # list 
                         chrlen.file = 'path to file with chromosome lengths;
                         set to NULL to not check',
                         bin.len = 'contact map resolution',
                         invalidij.action = 'specify what to do with unwanted
                         contacts; NA - set their values to NA, "drop" - remove
                         them from df, "none" - keep them as it is',
                         species.id = 'ath | dme | osa | human', 
                         categoriseValues = "before" # or "after" filtering
                         ){

  print(metric.dir)
  
  # Get df of upper triangle contacts
  if( metric%in%c("Cs.raw", "Cs.norm", "Cp") ){
    
    load(file=paste0(metric.dir, "/", species.id, "_", chr, "_allcontacts.RData"))
    # MELT.MX$upper.tri contains only contacts present in at least 1 tissue (with Cp value)
    MELT.MX$rest <- NULL
    ct.v <- colnames(MELT.MX$upper.tri[,-(1:2), drop=F])
    
    if(metric=="Cp"){
      
      value <- apply(X=MELT.MX$upper.tri[,-(1:2)], MARGIN=1, 
                     FUN=function(rw)sum(rw!=0))
      
      if(ct%in%ct.v){
        
        value[ MELT.MX$upper.tri[[ct]]==0 ] <- 0
        print("getContactDF(): Cp per cell/tissue...", quote=FALSE)
        
      } else if(ct=="All"){
        print("getContactDF(): Cp all cell/tissue...", quote=FALSE)
      } else {
        stop("getContactDF(): Invalid ct for Cp-based metric.")
      }
      
      MELT.MX$upper.tri <- cbind.data.frame(MELT.MX$upper.tri[,c("i","j")], 
                                            value=value)
      
      print("getContactDF(): Cp value obtained.", quote=F)
      
      rm(value)
      
    } else {
      
      if(!ct%in%ct.v){
        stop("getContactDF(): Invalid ct for Cs-based metric.")
      }
      MELT.MX$upper.tri <- MELT.MX$upper.tri[,c("i", "j", ct)]
      colnames(MELT.MX$upper.tri) <- colnames(MELT.MX$upper.tri.nocontact)
      
      print("getContactDF(): Cs.norm value obtained.", quote=F)
      
    }
    
    df <- rbind(MELT.MX$upper.tri, MELT.MX$upper.tri.nocontact)
    
    if( any(!is.finite(df$value)) | any(df$value < 0) ){
      warning("getContactDF(): Non-finite/Negative ", metric, " value.", quote=F)
    }
    
    rm(MELT.MX)
    
  } else if( grepl(x=metric, pattern="CII.cont.|CII.disc.") ){
    
    m.v <- strsplit(x=metric, split=".", fixed=TRUE)[[1]]
    if( length(m.v)!=4 ){
      stop("getContactDF(): Invalid metric argument for c|| values.")
    }
    names(m.v) <- c("metric", "value", "type", "cutoff")
    
    type.id <- m.v[["type"]]
    value <- "C||"
    if(type.id=="G"){
      
      type.id <- "kmer"
      value <- "Gfree"
      
    }
      
    #load(file=paste0(metric.dir, "/", chr, "_", type.id, "_", gcb, 
    #                 "_grouped_cutoff", m.v["cutoff"], ".RData"))
    load(file=paste0(metric.dir, "/", chr, "_", type.id, "_", gcb, ".RData"))
   
    # Determine which column of values to take
    df <- data.frame(CII.MX[,c("i", "j")], value=CII.MX[,value], stringsAsFactors=FALSE)
    rm(CII.MX)
    
    if( m.v[["value"]]=="disc" & categoriseValues == "before" ){
      
      if(m.v[["type"]]=="G"){
        # Negate G because G is negatively correlated with CII
        df$value <- -(df$value)
      } 
      
      df$value <- categoriseValues(val.v=df$value, cutoff=as.numeric(m.v[["cutoff"]]))
      
      print("getContactDF(): Categorising values before filtering.")
      
    }  
    
    print("getContactDF(): C|| value obtained.", quote=F)
      
    rm(type.id, value)
    
  } else if( grepl(x=metric, pattern="SIM.", fixed=TRUE) ){
    
    sim.mx <- data.matrix(fread(file=paste0(metric.dir, "/", gcb, "_", chr, "_sim.txt"), 
                                header=FALSE, data.table=FALSE, stringsAsFactors=FALSE))
    dimnames(sim.mx) <- NULL
    # Simulation will not include last bin of chr because it's shorter than resolution
    sim.mx <- cbind(rbind(sim.mx, NA), NA)
    len <- length(sim.mx[,1])
    if( !identical( dim(sim.mx), c(len,len) ) ){
      stop("getContactDF(): sim.mx not square matrix.")
    }
    dimnames(sim.mx) <- list(1:len, 1:len)
    df <- reshape2::melt(sim.mx, na.rm=FALSE)
    rm(sim.mx)
    colnames(df) <- c("i", "j", "value")
    # Take upper triangle
    df <- df[df$i<df$j,]
    
    print("getContactDF(): Csim value obtained.", quote=F)
    
  } else {
    stop(paste0("getContactDF(): Invalid ", metric, "!"))
  }
  
  gc()
  
  # Check that all columns should be numeric
  if( !is.numeric(c(df$i, df$j, df$value)) ){
    stop("getContactDF(): df has non-numeric column.")
  }
  
  # Check that all are upper triangle contacts
  if( any(df$i>=df$j) ){
    stop("getContactDF(): df not from upper triangle.")
  }
  
  if( !is.null(chrlen.file) ){
    
    # Check if number of contacts is correct based on number of bins per chr
    genome <- read.table(file=chrlen.file, stringsAsFactors=FALSE, header=TRUE)
                         #colClasses=c("character", "integer", "integer"))
    chr.bin.num <- ceiling(genome$length.bp[genome$chromosome==chr]/bin.len)
    #if(chr.bin.num!=genome$bins.40kb[genome$chromosome==chr]){
    #  stop("getContactDF(): chromosome length wrong.")
    #} 
    
    if( length(df$value)!=((chr.bin.num*chr.bin.num)-chr.bin.num)/2 ){
      stop("getContactDF(): Wrong number of upper triangle contacts.") 
    }
    
    rm(genome, chr.bin.num)
    
    print("getContactDF(): Number of contacts match expected value.")
    
  } 
 
  # Filter contacts based on specificed i and j bins of interest and gap range
  incl.TF <- filterContacts(ij.df=df[,c("i","j")], gap.range=gap.range,
                            incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y,  
                            mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y)
  df$include <- as.numeric(incl.TF)
  if( any(is.na(df$include)) ){
    stop("getContactDF(): NA in df$include.")
  }
  
  # Check for NA values when metric is Cs or Cp
  if( any(is.na(df$value)) & !grepl(x=metric, pattern="CII.|SIM.") ){
    warning("getContactDF(): NA values in df.")
  }
  
  # Check if all non-NA values are >=0 except for CII
  if( !grepl(x=metric, pattern="CII.") & any(df$value < 0 & !is.na(df$value)) ){
    warning("getContactDF(): Negative contact value.")
  }
  
  # Deal with unwanted/invalid contacts 
  if( is.na(invalidij.action) ){
    
    df$value[df$include==0] <- NA
    print("getContactDF(): Unwanted contacts' values set to NA.")
    
  } else if( !is.na(invalidij.action) & invalidij.action=="drop"){
    
    df <- df[df$include==1,]
    print("getContactDF(): Unwanted contacts removed.")
    
  } else if( !is.na(invalidij.action) & invalidij.action=="none"){
    print("getContactDF(): Unwanted contacts untouched.")
  } else {
    stop("getContactDF(): Invalid invalidij.action argument.")
  }
  
  # Categorise values after filtering
  
  if( grepl(x=metric, pattern="CII.disc.") ){
    
    if( m.v[["value"]]=="disc" & categoriseValues == "after" ){
      
      if(m.v[["type"]]=="G"){
        # Negate G because G is negatively correlated with CII
        df$value <- -(df$value)
      } 
      
      df$value <- categoriseValues(val.v=df$value, cutoff=as.numeric(m.v[["cutoff"]]))
      
      print("getContactDF(): Categorising values after filtering.")
      
    }
    
  }
  
  print(paste0("getContactDF(): ", metric, " values retrieved!"), quote=FALSE)
  
  return(df)
  
}
################################################################################
getContactDF <- cmpfun(getContactDF, options=list(suppressUndefined=TRUE))
################################################################################
# rm(list=ls()); gc()
