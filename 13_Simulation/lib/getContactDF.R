################################################################################
# Get dataframe of upper triangle contacts (i-j column) including a column for 
# values based on metric and a column specifying whether the contact will be 
# included or not based on specified bins to be included/masked
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(data.table)
# library(reshape)
### FUNCTION ###################################################################
getContactDF <- function(metric.dir=metric.dir, metric="Cp", gcb = "min2Mb", 
                         chr = "chr1", ct = "FC", gap.range=gap.range,
                         incl.bin.x = NULL, incl.bin.y = NULL, # list, 
                         mask.bin.x = NULL, mask.bin.y = NULL  # list 
                         ){
  
  p.lst <- list()
  if( metric%in%c("Cs.raw", "Cs.norm", "Cp") ){
    
    load(file=paste0(metric.dir, "/human_", chr, "_allcontacts.RData"))
    MELT.MX$rest <- NULL
    
    if(metric=="Cp"){
      value <- apply(X=MELT.MX$upper.tri[,-(1:2)], MARGIN=1, 
                     FUN=function(rw)sum(rw!=0))
      if( ct%in%colnames(MELT.MX$upper.tri[,-(1:2)]) ){
        drop.TF <- MELT.MX$upper.tri[[ct]]==0
        value[drop.TF] <- 0; rm(drop.TF)
        print("Cp per cell/tissue...", quote=FALSE)
      } else if(ct=="hg19"){
        print("Cp all cell/tissue...", quote=FALSE)
      } else {
        stop("Invalid ct.")
      }
      MELT.MX$upper.tri <- cbind.data.frame(MELT.MX$upper.tri[,c("i","j")], 
                                            value=value); rm(value)
    } else{
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
    diag.val <- unique(diag(sim.mx))
    if( length(diag.val)!=1 ){ stop("Diag value in sim.mx not unique.") }
    # Simulation will not include last bin of chr because it's shorter than resolution
    sim.mx <- cbind(rbind(sim.mx, NA), NA)
    #if( !isSymmetric(sim.mx) ){ stop("sim.mx not symmetric.") }
    df <- reshape::melt(sim.mx); rm(sim.mx)
    colnames(df) <- c("i", "j", "value")
    # Take upper triangle
    df <- df[df$i<df$j,]
    df$value <- df$value/diag.val; rm(diag.val)
  
  } else {
    stop(paste0(metric, ": invalid!"))
  }
  
  gc()
  
  if( any(df$i>=df$j) ){ stop("df not from upper triangle.") }
  
  # Filter contacts based on specificed i and j bins; marked with -Inf value
  incl.TF <- filterContacts(ij.df=df[,c("i","j")], gap.range=gap.range,
                            incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y,  
                            mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y)
  df$include <- as.numeric(incl.TF)
  
  print(paste0(metric, " done!"), quote=FALSE)
  
  return(df)
  
}
################################################################################
# rm(list=ls()); gc()
