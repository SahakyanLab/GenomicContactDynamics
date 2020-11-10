################################################################################
# UniqueIDs of longest transcripts per contact persistence
# Unique ID is equivalent to row number in annotation table containing all 
# refseq categories
# For gene length analyses, only coding transcripts (NM_) will be considered
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
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    annofile.dir = "/t1-data/user/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
genelist.dir = paste0(wk.dir, "/out_anno_union")
out.dir = paste0(wk.dir, "/out_LTrUniqueIDPerCp")
### OTHER SETTINGS #############################################################
# gcb 
gcb.v = c("min2Mb", "min05Mb")
refseq.v = c("ALL", "NM", "NR")
# Annotation file prefix
anno.nme = "hg19anno"
# Text file lines
nCPU = 4
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools) # isplitVector
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# This should be the table loaded regardless of refseq category because the 
# uniqueID for each transcript is based on its row number in this table
annotable <- fread(file=paste0(annofile.dir, "/", anno.nme, "_ALL"), 
                   header=TRUE, data.table=FALSE, stringsAsFactors=FALSE) 

annotable$txSize <- annotable$txEnd-annotable$txStart
#if(!all(annotable$txSize)>0){
#  stop("Negative transcript lengths.")
#}

# Confirm that uniqueID column in annotation table correspond to row numbers
# so we can treat uniqueID as row numbers.
#if( !identical( as.integer(rownames(annotable)), annotable$uniqueID) ){
#  stop("UniqueID not increasing (not equivalent to row numbers).")
#}

for(gcb in gcb.v){
  for(refseq in refseq.v){
    uniqueIDperCpStr <- readLines(con=paste0(genelist.dir, "/", gcb, "_", refseq, "_uniqueID"))
    # Exclude entries for "HiC_all"
    HiC_all.ind <- grep(x=uniqueIDperCpStr, pattern="HiC", fixed=TRUE)
    uniqueIDperCpStr <- uniqueIDperCpStr[-(c(HiC_all.ind, HiC_all.ind+1L))]
    
    
    header.ind <- grep(x=uniqueIDperCpStr, pattern=">", fixed=TRUE)
    uniqueIDs <- uniqueIDperCpStr[header.ind+1]
    uniqueIDs.len <- length(uniqueIDs)
    
    LtrUniqueIDstring <- foreach(itr=isplitVector(x=1:uniqueIDs.len , chunks=nCPU),
                                 .inorder=TRUE, .combine="c",
                                 .export=c("annotable", "uniqueIDs"), 
                                 .noexport=ls()[!ls()%in%c("annotable", "uniqueIDs")]
    ) %op% {
      
      chunk <- sapply(X=itr, FUN=function(i){
        genestring <- uniqueIDs[i]
        uniqueID <- as.numeric( strsplit(x=genestring, split=";")[[1]] )
        
        # Subset annotation table with desired uniqueIDs
        annotable.sub <- annotable[uniqueID, c("uniqueID", "txSize", "name2", "name")]
        rm("genestring", "uniqueID", "annotable"); gc()
        
        # UniqueID of longest transcripts per cp/strength
        LtrPerHugo <- by(data=annotable.sub[, c("uniqueID", "txSize", "name")], 
                         INDICES=annotable.sub$name2, 
                         FUN=function(df){
                           # Reduce name (refseq ID) to NM and NR only so 
                           # ordering by name then unique ID is possible
                           df$name <- gsub(x=df$name, replacement="",
                                           pattern="[_].*$")
                           # Order (increasing in number and alphabet) based on 
                           # name then unique ID. This is for picking 
                           # from top of equal length transcripts and prefentially
                           # picking NM over NR in case both are present. 
                           df <- df[order(df$name, df$uniqueID, decreasing=FALSE),]
                           df[df$txSize==max(df$txSize, na.rm=TRUE),
                              "uniqueID"][1]
                         })
        # String of uniqueIDs of longest transcripts per cp
        paste(LtrPerHugo, collapse=";")
      })
      
    }
    
    uniqueIDperCpStr[header.ind+1L] <- LtrUniqueIDstring
    
    write(x=uniqueIDperCpStr, file=paste0(out.dir, "/", gcb, "_", refseq, 
                                          "_LtrUniqueIDstringPerCpStr"))
    
    rm(uniqueIDs, uniqueIDperCpStr, LtrUniqueIDstring); gc()
  }
}

# rm(list=ls())
