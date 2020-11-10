################################################################################
# Make a file of Cs (per tissue) and Cp values per minimum distance between
# contacting regions (min2Mb/min05Mb)
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
    persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    persist.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", 
          quote=FALSE)
  }
}
out.dir = paste0(wk.dir, "/out_CsCpValues")
### OTHER SETTINGS #############################################################
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
gcb.v = c("min2Mb", "min05Mb")
chr.v = paste("chr", c(1:22, "X"), sep="")
nCPU = 2L
out.name = "CsCpValues"
### FUNCTION ###################################################################
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools) #isplitVector
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)
celltiss.v <- c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
                "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC", "cp")
celltiss.v.len <- length(celltiss.v)

toExport <- c("persist.dir", "out.dir", "out.name", "chr.v", "chr.v.len", 
              "celltiss.v", "celltiss.v.len", "lst")

#### PARALLEL EXECUTION #########

foreach(gcb=gcb.v, .inorder=FALSE, .export=toExport, 
        .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in 1:chr.v.len){
    
    chr <- chr.v[i]
    load( file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData") )
    PERSIST.MX$hits[["cp"]] <- PERSIST.MX[["ntis"]]
    print(chr)
    
    if(i!=1){
      
      lst.temp <- sapply(X=celltiss.v, simplify=FALSE, FUN=function(id.nme){
        unique(PERSIST.MX$hits[[id.nme]])
      })
      lst <- sapply(X=celltiss.v, simplify=FALSE, FUN=function(id.nme){
        union(lst[[id.nme]], lst.temp[[id.nme]])
      })
      rm(lst.temp); gc()
      
    } else {

      lst <- sapply(X=celltiss.v, simplify=FALSE, FUN=function(id.nme){
        unique(PERSIST.MX$hits[[id.nme]])
      })
      print("lst initiated.")
      names(lst) <- celltiss.v
    }
  
  } # chr.v for loop end
  
  out <- sapply(X=names(lst), simplify=FALSE, FUN=function(id.nme){
    vals <- sort(unique(lst[[id.nme]]), decreasing=FALSE)
    return( c(id.nme, paste(vals, collapse=";")) )
  })
  
  write.table(x=do.call("rbind", out), 
              file=paste0(out.dir, "/", gcb, "_", out.name), 
              col.names=c("id", "values"),
              row.names=FALSE, sep="\t", quote=FALSE)
} 

### END OF PARALLEL EXECUTION ###

# rm(list=ls())

  



