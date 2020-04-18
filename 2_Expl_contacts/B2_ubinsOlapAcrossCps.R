################################################################################
# Determine and plot the percent overlap of unique contacting bins across Cps.
# Unique bins from all chromosomes were combined.   
# Overlaps tested: HiCAll_1, 1_2...1_21. Values are percent overlap. For example,
# 1_21==10% for FC cell line means that 10% of all Cp=21 unique bins are also in
# Cp=1 unique bins.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = persist.dir
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb.v = c("min2Mb", "min05Mb")
chr.v = paste0("chr", c(22:1, "X"), sep="")
nCPU = 10L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
ct.v <- c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
          "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
ct.v.len <- length(ct.v)
chr.v.len <- length(chr.v)

query.v <- as.character(c(1:21))
query.v.len <- length(query.v)
ref.v <- c("HiCAll", rep(x=1, times=20))

chrLen.df <- fread(file=chrLenfile, colClasses=list(character=1, integer=2), 
                   stringsAsFactors=FALSE, data.table=FALSE)

for(gcb in gcb.v){
  
  toExport <- c("ct.v", "chr.v", "persist.dir", "gcb", "query.v", 
                "query.v.len", "ref.v")
  
  #### PARALLEL EXECUTION #########
  
  OLAPCP.MX <- foreach(itr=isplitVector(1:ct.v.len, chunks=nCPU), 
                       .inorder=TRUE, .combine="rbind",
                       .export=toExport, 
                       .noexport=ls()[!ls()%in%toExport]
                       
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      ct <- ct.v[i]
      
      # Collect bins per chromosome separating per Cp
      binPCP.lst <- vector(mode="list", length=21L)
      names(binPCP.lst) <- 1:21
      
      for(c in 1:chr.v.len){
        
        chr <- chr.v[c]
        chr.num <- gsub(x=chr, pattern="chr", replacement="", fixed=TRUE)
        
        # Load PERSIST.MX directory
        load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
        log <- PERSIST.MX$hits[,ct]!=0L
        
        # Per Cp list of bins written as Cp_bin, e.g. "22_503"
        temp.lst <- by(data=PERSIST.MX$hits[log,c("i", "j")],
                       INDICES=PERSIST.MX$ntis[log], 
                       FUN=function(x){
                         x <- unique(unlist(x))
                         if(any(is.na(x))){
                           stop("Checkpoint 1.")
                         }
                         return( paste(chr.num, x, sep="_") )
                       })
        
        rm(log)
        
        cp.v <- names(temp.lst)
        binPCP.lst <- lapply(X=cp.v, FUN=function(cp){
          unique( c(binPCP.lst[[cp]], temp.lst[[cp]]) )
        })
        names(binPCP.lst) <- cp.v
        
        rm(PERSIST.MX, temp.lst, cp.v, chr.num); gc()
        
        print(paste0(chr, " data added!"), quote=FALSE)
        
      } # chr.v.len for loop end
      
      percOlap.v <- rep(NA, times=query.v.len)
      for(k in 1:query.v.len){
        
        query <- query.v[k]
        ref <- ref.v[k]
        
        if(ref!="HiCAll"){
          olap <- intersect( binPCP.lst[[ref]], binPCP.lst[[query]] )
        } else if(ref=="HiCAll"){
          olap <- intersect( unique(unlist(binPCP.lst)), binPCP.lst[[query]] )
        } else {
          stop("Invalid ref.")
        }
        percOlap.v[k] <- 100*( length(olap)/length(binPCP.lst[[query]]) )
        
      } # olap.v for loop end
      
      rm(binPCP.lst); gc()
      
      print(paste0(ct, " done!"), quote=FALSE)
      
      return(percOlap.v)
      
    })
    
    return(do.call("rbind",  chunk))
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  dimnames(OLAPCP.MX) <- list(ct.v, paste(ref.v, query.v, sep="_"))
  save(OLAPCP.MX, file=paste0(out.dir, "/chrALL_", gcb, "_ubinsOlapCp.RData"))
  
} # gcb.v for loop end

# rm(list=ls())



 

