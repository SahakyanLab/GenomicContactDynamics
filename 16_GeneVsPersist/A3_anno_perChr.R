################################################################################
# Per chr, identify genes co-localising with contacts per Cp category 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/16_GeneVsPersist"
    persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
    persist.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
feat.bin.dir = paste0(wk.dir, "/out_mapToPersistBins_anno")
out.dir = paste0(wk.dir, "/out_anno_perChr")
### OTHER SETTINGS #############################################################
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
gcb = "min2Mb" # "min2Mb" | "min05Mb"
refseq.v = c("LTr_ALL", "LTr_NM", "LTr_NR")
chr.v = paste("chr", c(22:1, "X"), sep="")
nCPU = 3L
cutoffStr = 5L
### FUNCTION ###################################################################
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools) #isplitVector
source(paste0(lib, "/UTL_doPar.R"))

# Don't use comma as splitChar
uniqStringVec <- function(vec=FEATURE.BIN.MX[FEATURE.BIN.MX$bin%in%bins.uniq, "name2"],
                          splitChar=";"){
  x <- unlist(x=strsplit(x=vec, split=splitChar), use.names=FALSE)
  return( paste(unique(x), collapse=splitChar) )
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)
celltiss.v <- c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
              "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
toExport <- c("feat.bin.dir", "chr.v", "celltiss.v", "out.dir", 
              "cutoffStr", "gcb")

for(refseq in refseq.v){
  
  toExport <- unique(c(toExport, "refseq"))
  
  #### FOREACH EXECUTION #########
  
  foreach( itr=isplitVector(1:chr.v.len, chunks=nCPU), .inorder=FALSE,
           .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    
    for(i in itr){
     
      chr <- chr.v[i]
      
      # Load FEATURE.BIN.MX
      load( file=paste0(feat.bin.dir, "/", 
                        dir(path=feat.bin.dir, 
                            pattern=paste0(chr, "_", gcb, "_", refseq), full.names=FALSE)) )
      FEATURE.BIN.MX <- FEATURE.BIN.MX[FEATURE.BIN.MX$countPerBin>0,]
      
      #---------------CONTACTPERSISTENCE------------
      
      # Load PERSIST.MX, contacts are unique
      load( file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData") )
      ntis.v <- 1:21
      
      for(ntis in ntis.v){
       
        # subset contacts by ntis
        binsperNtis <- PERSIST.MX$hits[PERSIST.MX$ntis==ntis, c("i", "j")]
        bins.uniq <- unique( c(unique(binsperNtis[,"i"]),
                               unique(binsperNtis[,"j"])) )
        
        rm("binsperNtis"); gc()
        
        header <- paste0(">all_genes_cp_", ntis, "_end")
        sapply(X=c("name2", "uniqueID"), FUN=function(x){
          
          x.uniq <- uniqStringVec(vec=FEATURE.BIN.MX[FEATURE.BIN.MX$bin%in%bins.uniq, x],
                                  splitChar=";")
          filename <- paste0(out.dir, "/", chr, "_", gcb, "_", refseq, 
                             "_", x)
          write(x=header, file=filename, append=TRUE)
          write(x=x.uniq, file=filename, append=TRUE)
        }) # sapply
        
      } # ntis end for loop
      
      #---------------CONTACTHICSTRENGTH------------
      
     allTissPerChr <- sapply(X=celltiss.v, simplify=FALSE,
                             FUN=function(tiss){
        
        grthanCutoff <- paste0("grthan", cutoffStr)
        count.v <- c("HiC_all", 1:cutoffStr, grthanCutoff)
        
        allStrPerTiss <- sapply(X=count.v, simplify=FALSE, 
                                FUN=function(strength){
          
          # subset contacts by tissue and strength
          if(strength=="HiC_all"){
            binsperTissStrength <- PERSIST.MX$hits[PERSIST.MX$hits[[tiss]]>0, 
                                                   c("i", "j")]
          } else if(strength==grthanCutoff){
            binsperTissStrength <- PERSIST.MX$hits[PERSIST.MX$hits[[tiss]]>cutoffStr, 
                                                   c("i", "j")]
          } else {
            binsperTissStrength <- PERSIST.MX$hits[PERSIST.MX$hits[[tiss]]==strength, 
                                                   c("i", "j")]
          }
          
          bins.uniq <- unique( c(unique(binsperTissStrength[,"i"]),
                                 unique(binsperTissStrength[,"j"])) )
          
          rm("binsperTissStrength"); gc()
          
          test <- (FEATURE.BIN.MX$bin)%in%bins.uniq
          c(paste0(">all_genes_", tiss, "_count_", strength, "_end"), 
            uniqStringVec(vec=FEATURE.BIN.MX[test, "name2"],
                          splitChar=";"),
            uniqStringVec(vec=FEATURE.BIN.MX[test, "uniqueID"],
                          splitChar=";"))
         
        })  # sapply cutoff
        
        unlist(x=allStrPerTiss, use.names=FALSE)
        
      }) # sapply tissue
      
      all.vec <- unlist(allTissPerChr, use.names=FALSE)
      all.vec.len <- length(all.vec)
      num <- length(celltiss.v)*(cutoffStr+2L)*3
    
      if(all.vec.len==num){
       
        write(x=all.vec[-(seq(from=3, to=all.vec.len, by=3))], 
              file=paste0(out.dir, "/", chr, "_", gcb, "_", refseq,
                          "_name2"),
              append=TRUE)
        
        write(x=all.vec[-(seq(from=2, to=(all.vec.len-1L), by=3))], 
              file=paste0(out.dir, "/", chr, "_", gcb, "_", refseq,
                          "_uniqueID"),
              append=TRUE)
      }
      
      rm("PERSIST.MX", "FEATURE.BIN.MX", "allTissPerChr"); gc()  
      
      print(paste0(chr, " done!"))
      
    } # itr for loop end
    
  } # foreach loop end
  
  ### END OF FOREACH EXECUTION ###
  
} # refseq.v for loop end

# rm(list=ls()); gc()

  



