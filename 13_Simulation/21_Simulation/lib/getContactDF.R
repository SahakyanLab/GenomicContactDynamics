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
# source(paste0(wk.dir, "/lib/filterContacts.R"))
### FUNCTION ###################################################################
getContactDF <- function(metric.dir=metric.dir, metric="Cp", gcb = "min2Mb", 
                         chr = "chr1", ct = "FC", gap.range=gap.range,
                         incl.bin.x = NULL, incl.bin.y = NULL, # list, 
                         mask.bin.x = NULL, mask.bin.y = NULL,  # list 
                         chrlen.file = 'path to file with chromosome lengths;
                         set to NULL to not check',
                         bin.len = 'contact map resolution',
                         invalidij.action = 'specify what to do with unwanted 
                         contacts; NA - set their values to NA, "drop" - remove
                         them from df, "none" - keep them as it is'
                         ){
  
  p.lst <- list()
  
  # Get df of upper triangle contacts
  if( metric%in%c("Cs.raw", "Cs.norm", "Cp") ){
    
    load(file=paste0(metric.dir, "/human_", chr, "_allcontacts.RData"))
    # MELT.MX$upper.tri contains only contacts present in at least 1 tissue (with Cp value)
    MELT.MX$rest <- NULL
    ct.v <- colnames(MELT.MX$upper.tri[,-(1:2)])
    
    if(metric=="Cp"){
      
      value <- apply(X=MELT.MX$upper.tri[,-(1:2)], MARGIN=1, 
                     FUN=function(rw)sum(rw!=0))
      
      if(ct%in%ct.v){
        
        value[ MELT.MX$upper.tri[[ct]]==0 ] <- 0
        print("getContactDF(): Cp per cell/tissue...", quote=FALSE)
        
      } else if(ct=="hg19"){
        print("getContactDF(): Cp all cell/tissue...", quote=FALSE)
      } else {
        stop("getContactDF(): Invalid ct for Cp-based metric.")
      }
      
      MELT.MX$upper.tri <- cbind.data.frame(MELT.MX$upper.tri[,c("i","j")], 
                                            value=value)
      rm(value)
      
    } else {
      
      if(!ct%in%ct.v){
        stop("getContactDF(): Invalid ct for Cs-based metric.")
      }
      MELT.MX$upper.tri <- MELT.MX$upper.tri[,c("i", "j", ct)]
      colnames(MELT.MX$upper.tri) <- colnames(MELT.MX$upper.tri.nocontact)
      
    }
    
    MELT.MX$upper.tri <- rbind(MELT.MX$upper.tri, MELT.MX$upper.tri.nocontact)
    # Take upper triangle
    df <- MELT.MX$upper.tri
  
    rm(MELT.MX)
    
  } else if( grepl(x=metric, pattern="CII.cont.|CII.disc.") ){
    
    m.v <- strsplit(x=metric, split=".", fixed=TRUE)[[1]]
    names(m.v) <- c("metric", "value", "type", "cutoff")
    load(file=paste0(metric.dir, "/", chr, "_", m.v["type"], "_", gcb, 
                     "_grouped_cutoff", m.v["cutoff"], ".RData"))
    value <- ifelse(m.v["value"]=="disc", "group", "C||")
    # Take upper triangle
    df <- cbind.data.frame(CII.MX[,c("i", "j")], value=CII.MX[,value])
    rm(CII.MX, m.v, value)
    
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
    genome <- read.table(file=chrlen.file, stringsAsFactors=FALSE, header=TRUE,
                         colClasses=c("character", "integer", "integer"))
    chr.len <- ceiling(genome$length.bp[genome$chromosome==chr]/bin.len)
    if(chr.len!=genome$bins.40kb[genome$chromosome==chr]){
      stop("getContactDF(): chromosome length wrong.")
    } 
    
    if( length(df$value)!=((chr.len*chr.len)-chr.len)/2 ){
      stop("getContactDF(): Wrong number of upper triangle contacts.") 
    }
    
    rm(genome, chr.len)
    
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
    stop("getContactDF(): NA values in df.")
  }
  
  # Check if all non-NA values are >=0 except for CII
  if( !grepl(x=metric, pattern="CII.") & any(df$value<0 & !is.na(df$value)) ){
    stop("getContactDF(): Negative contact value.")
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
  
  print(paste0("getContactDF(): ", metric, " values retrieved!"), quote=FALSE)
  
  return(df)
  
}
################################################################################
getContactDF <- cmpfun(getContactDF, options=list(suppressUndefined=TRUE))
################################################################################
# rm(list=ls()); gc()
