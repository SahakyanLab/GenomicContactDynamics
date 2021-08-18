################################################################################
# Find regions forming contacts of desired Cs to Cp correlation
# (e.g. directly/inversely proportional) 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/4_CsVsCp"
    data.dir = "/Users/ltamon/Database"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/13_CsVsCp"
    data.dir = "/t1-data/user/ltamon/Database"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_findRegion")
### OTHER SETTINGS #############################################################
#ct.v = c("FC", "LC", "ESC", "Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", 
#         "Li", "SB", "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC")
ct.v = "FC"
gcb.v = c("min2Mb")
chr.v = "chr17" #paste("chr", c(1:22, "X"), sep="")
# Range of Cs, [x,y]
Cs = c(1,3)
# Range of Cp, [x,y]
Cp = c(1,5)
out.name = "1and2CsLowCp"
scaled = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
### FUNCTION ###################################################################
getContactBasedOnCsCpRange <- function(
  out.dir = "/dir",
  out.name = "lowCshighCp",
  PERSIST.MX = PERSIST.MX,
  gcb = "min2Mb",
  chr = "chr21",
  ct = "FC", 
  Cs = c(1,4),
  Cp = c(19,21)
){
  
  log.ct <- PERSIST.MX$hits[[ct]]!=0
  ij.mx <- cbind(PERSIST.MX$hits[log.ct,c("i", "j", ct)],
                 Cp=PERSIST.MX$ntis[log.ct])
  rm(PERSIST.MX, log.ct); gc()
  
  # Identify contacts satisfying Cs and Cp desired ranges
  ij.mx <- ij.mx[ ij.mx[,ct]%in%(Cs[1]:Cs[2]) & ij.mx[,"Cp"]%in%(Cp[1]:Cp[2]), ]
  # Make sure i is in increasing order
  ij.mx <- ij.mx[order(ij.mx[,"i"], decreasing=FALSE),]
  
  write.table(ij.mx, file=paste0(out.dir, "/", chr, "_", gcb, "_", ct, "_", out.name),
              quote=FALSE, sep="\t", row.names=FALSE)
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Old lengths
len.ct<- length(ct.v)
len.gcb <- length(gcb.v)
len.chr <- length(chr.v)

# New vectors
ct.v <- rep(ct.v, each=len.gcb*len.chr)
gcb.v <- rep(gcb.v, each=len.ct*len.chr)
chr.v <- rep(chr.v, each=len.ct*len.gcb)

# New lengths
len.ct <- length(ct.v)
len.gcb <- length(gcb.v)
len.chr <- length(chr.v)

if(unique(len.ct, len.gcb, len.chr)==len.ct){
  
  affix <- ifelse(scaled==TRUE, "_scaled", "")
  
  for(i in 1:len.ct){
    
    ct <- ct.v[i]
    gcb <- gcb.v[i]
    chr <- chr.v[i]
    
    # Load PERSIST.MX (original/scaled)
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, affix, ".RData"))
    
    getContactBasedOnCsCpRange(
      out.dir=out.dir,
      out.name=paste0(out.name, affix),
      PERSIST.MX=PERSIST.MX,
      gcb=gcb,
      chr=chr,
      ct=ct,
      Cs=Cs,
      Cp=Cp
    )
    
    print(paste0(gcb, "_", chr, "_", ct, " done!"), quote=FALSE)
    
  } # len.ct for loop end
  
}

# rm(list=ls())