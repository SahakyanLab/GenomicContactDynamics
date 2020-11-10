################################################################################
# Combine TRCOUNT.MX of all chr.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    objective.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    objective.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
output.dir = paste0(objective.dir, "/out_")
### OTHER SETTINGS #############################################################
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
gcb = "min2Mb" #min05Mb
id.v = c("cp", "Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
         "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
refseq = "ALL" 
chr.v = paste("chr", c(1:22, "X"), sep="")
nCPU = 23L
Lref.v = "ave" #c("orig", "ave", "m2sd")
### FUNCTION ###################################################################
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(itertools) #isplitVector
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/loadRData.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)
print(gcb)

for(l in Lref.v){
  
  out.name <- paste0("ContactCountPerTr_", l)
  toExport <- c("output.dir", "chr.v", "out.name")
  
  for(id in id.v){
    
    toExport <- unique(c(toExport, "id"))
    print(id)
    #### PARALLEL EXECUTION #########
    
    TRCOUNT.MX <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                          .inorder=TRUE, .combine="rbind",
                          .export=toExport, .noexport=ls()[!ls()%in%toExport]
    ) %op% {
      
      chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
        chr <- chr.v[i]
        affix <- paste0(chr, "_", gcb, "_", refseq, "_", id)
        mx <- loadRData(file=paste0(output.dir, "/", affix, "_", out.name, ".RData"))
        print(chr)
        mx
      })
      
      do.call("rbind", chunk)
      
    }
    ### END OF PARALLEL EXECUTION ###
    
    affix <- paste0("chrALL_", gcb, "_", refseq, "_", id)
    save(TRCOUNT.MX, file=paste0(output.dir, "/", affix, "_", out.name, ".RData"))
    
  } # id.v for loop end
  
  print(out.name)
  
} # Lref.v for loop end

# rm(list=ls())


  



