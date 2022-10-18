################################################################################
# Generate .RData per tissue with contact, Cf/Cs and Cp information
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/4_CsVsCp")
out.dir = paste0(wk.dir, "/out_CsVsCpobj")
### OTHER SETTINGS #############################################################
ct.v = "REPLACE" #sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
      #        "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "LC", "FC"))
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
nCPU = 1L #~15G
bin.len = 40000

persist.id = "rawCs"
if(persist.id == "rawCs"){
  persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/")
} else if(persist.id == "HiCNormCs"){
  persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm")
} else{
  stop("Invalid persist.id. Either 'rawCs' or 'HiCNormCs'.")
}
print(paste0("persist.dir: ", persist.dir))
out.dir = paste0(out.dir, "/", persist.id)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)
toExport <- c("chr.v", "persist.dir", "gcb")

for(ct in ct.v){
  
  id <- paste0("chrALL_", gcb, "_", bin.len, "bpbin_", ct)
  
  #### PARALLEL EXECUTION #########
  CPCS.MX <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                     .inorder=T, .combine="rbind",
                     .export=toExport, .noexport=ls()[!ls()%in%toExport]
                     
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=F, FUN=function(i){
      
      chr <- chr.v[i]
      load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      rownames(PERSIST.MX$hits) <- NULL
      
      gaps.bin <- PERSIST.MX$hits$j - PERSIST.MX$hits$i - 1
      ct.TF <- PERSIST.MX$hits[,ct] != 0
      
      return(
        cbind(Cs=PERSIST.MX$hits[ct.TF,ct], 
              Cp=PERSIST.MX$ntis[ct.TF],
              gap.jminusiminus1.bin=gaps.bin[ct.TF])
      )
      
    })
    
    print(paste0(ct, " done!"), quote=F)
    return( do.call("rbind", chunk) ) 
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  save(CPCS.MX, file=paste0(out.dir, "/", id, "_CsVsCp.RData"))
  
  rm(CPCS.MX, id)
  gc()
  
}

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()