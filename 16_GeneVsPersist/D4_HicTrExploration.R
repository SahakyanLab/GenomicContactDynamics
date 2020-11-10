################################################################################
# Per transcript, count the number of its overlapping contacts having Cs/Cp 
# values per tissue
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
annofile.dir = paste0(data.dir, "/ucsc_tables/hsa_geneAnno")
feat.bin.dir = paste0(wk.dir, "/out_eq_mapToPersistBins_anno/FEATUREBINMX")
out.dir = paste0(wk.dir, "/out_HiCTrExploration")
### OTHER SETTINGS #############################################################
gcb = "min2Mb" # "min2Mb" | "min05Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
id.v = c("cp", "Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
         "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
refseq = "ALL" # "ALL" | "NM" | "NR")
Lref.v = c("orig", "ave", "m2sd")
# Max is total number of contacts
nCPU = 4L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools) #isplitVector
library(compiler)
source(paste0(lib, "/UTL_doPar.R"))
### FUNCTION ###################################################################
HicTrExploration <- function(
feat.bin.dir = paste0(wk.dir, "/out_eq_mapToPersistBins_anno/FEATUREBINMX_orig"),
out.dir = paste0(wk.dir, "/out_HiCTrExploration"),
gcb = "min2Mb",
id.v = c("cp", "Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
         "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC"),
refseq = "ALL",
chr = "chr21",
# Max is total number of contacts
nCPU = 1L,
out.name = "ContactCountPerTr_orig"
){

  # Load FEATURE.BIN.MX
  load( file=paste0(feat.bin.dir, "/", 
                    dir(path=feat.bin.dir, 
                        pattern=paste0(chr, "_", gcb, "_", refseq))) )
  feat.bin <- as.numeric(FEATURE.BIN.MX$bin)
  feat.bin.len <- length(feat.bin)
  
  bin.uniqueID <- strsplit(x=FEATURE.BIN.MX[,"uniqueID"], split=";")
  names(bin.uniqueID) <- feat.bin
  
  rm(FEATURE.BIN.MX); gc()
  
  bin.uniqueID <- foreach( bin=names(bin.uniqueID), 
                           .inorder=TRUE, .combine="rbind"
  ) %do% {
    return( cbind( bin=rep(as.numeric(bin)), 
                   uniqueID=as.numeric(bin.uniqueID[[bin]]) )
    )
  }
  
  bin.uniqueID <- by(data=bin.uniqueID[,"bin"],
                     INDICES=bin.uniqueID[,"uniqueID"],
                     FUN=function(x){unique(x)})
  uniqueID.v <- names(bin.uniqueID)
  uniqueID.v.len <- length(bin.uniqueID)
  
  # Load PERSIST.MX, contacts are unique
  load( file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData") )
  ij.mx <- cbind(PERSIST.MX$hits, cp=PERSIST.MX$ntis)
  
  # Remove contacts with no overlapping transcripts
  ij.mx <- ij.mx[ij.mx[,"i"]%in%feat.bin |
                   ij.mx[,"j"]%in%feat.bin, ]

  rm(PERSIST.MX); gc()
  
  # Possible values per id
  idVal.str <- fread(file=paste0(out.dir, "/", gcb, "_", "CsCpValues"),
                  header=TRUE, select=c("id", "values"), data.table=FALSE,
                  stringsAsFactors=FALSE)
  
  for(id in id.v){
   
    ij.mx.sub <- ij.mx[, c("i", "j", id)]
    
    #minIdVal <- 0
    
    #if(id=="cp"){
    #  maxIdVal <- 21L; minIdVal <- 1L
    #} else {
    #  maxIdVal <- maxval[maxval$celltiss==id, "maxCs"]; minIdVal <- 0
    #}
      
    toExport <- c("uniqueID.v", "bin.uniqueID", "ij.mx.sub", "id")
    
    #### PARALLEL EXECUTION #########
    
    TRCOUNT.MX <- foreach( itr=isplitVector(1:uniqueID.v.len, 
                                                chunks=nCPU), 
                               .inorder=TRUE, .combine="rbind",
                               .export=toExport, 
                               .noexport=ls()[!ls()%in%toExport]
    ) %op% {
      
      chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
        
        # Bins overlapping with transcript
        tr.bins <- bin.uniqueID[[uniqueID.v[i]]]
        
        # Count contact containing bins (per id)
        tr.id <- ij.mx.sub[ij.mx.sub[,"i"]%in%tr.bins |
                       ij.mx.sub[,"j"]%in%tr.bins, id]
        
        idVals <- strsplit(x=idVal.str[idVal.str$id==id, "values"], split=";")[[1]]
        idVals <- sort( as.numeric(idVals), decreasing=FALSE )
        tr.id.count <- table(tr.id)
        idVal.abse <- setdiff( idVals, as.numeric(names(tr.id.count)) )
        add.v <- rep(0, times=length(idVal.abse))
        names(add.v) <- idVal.abse
        tr.id.count <- c(tr.id.count, add.v)
        tr.id.count <- tr.id.count[ order( factor(names(tr.id.count), 
                                                     levels=idVals ) ) ]
        
      })
      
      do.call("rbind", chunk)
      
    }
    
    ### END OF PARALLEL EXECUTION ###
    
    rownames(TRCOUNT.MX) <- uniqueID.v
    affix <- paste0(chr, "_", gcb, "_", refseq, "_", id)
    save(TRCOUNT.MX, file=paste0(out.dir, "/", affix, "_", out.name, ".RData"))
    print(id)
    
    rm(TRCOUNT.MX, ij.mx.sub); gc()
    
  } # id.v for loop end 
} 
################################################################################
HicTrExploration <- cmpfun(HicTrExploration, options=list(suppressUndefined=TRUE))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)
chr.v <- rep(chr.v, times=length(Lref.v))
Lref.v <- rep(Lref.v, each=chr.v.len)
ext.len <- length(Lref.v)

for( refseq in refseq.v ){
  for(i in 1:ext.len){
    HicTrExploration(
      feat.bin.dir=paste0(feat.bin.dir, "_", Lref.v[i]),
      out.dir=out.dir,
      gcb=gcb,
      id.v=id.v,
      refseq=refseq,
      chr=chr.v[i], 
      nCPU=nCPU,
      out.name=paste0("ContactCountPerTr_", Lref.v[i])
    )
  }
} # Lref.v for loop end

# rm(list=ls())

  



