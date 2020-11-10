################################################################################
# 1-to-1 correspondence of gene with Cp (max Cp across transcript)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    objective.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    objective.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    annofile.dir = "/t1-data/user/ltamon/Database/ucsc_tables/hsa_geneAnno"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    objective.dir = "/home/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    annofile.dir = "/home/ltamon/Database/ucsc_tables/hsa_geneAnno"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
genelist.dir = paste0(objective.dir, "/out_anno_union")
output.dir = paste0(objective.dir, "/out_geneCentric")
### OTHER SETTINGS #############################################################
# gcb 
gcb = "min05Mb" #"min05Mb"
refseq.v = c("ALL", "NM", "NR")
# annotation file prefix
anno.nme = "hg19anno"
contactFeature = c("cp", "cs")
celltiss.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
               "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
nCPU=4L
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(itertools) # isplitVector
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# UniqueID corresponds to rows of this table of genes (name2)
annotable <- fread(file=paste0(annofile.dir, "/", anno.nme, "_ALL"), 
                   header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                   select=c("name2"))[[1]] 

identifier.v <- c()
if( !is.null(contactFeature) ){
  if("cp"%in%contactFeature){
    identifier.v <- "cp"
  }
  if("cs"%in%contactFeature){
    identifier.v <- c(identifier.v, celltiss.v)
  }
} else {
  stop("Invalid input for contactFeature (cp and/or cs).")
}

toExport <- c("annotable", "output.dir", "gcb")

for(refseq in refseq.v){
  
  # *_count_HiC_all entries already removed in this file
  LtrUniqueIDstringPerCpStr <- readLines(con=paste0(genelist.dir, "/", gcb, "_", 
                                                    refseq, "_LtrUniqueIDstringPerCpStr"))

  toExport <- unique(c(toExport, "refseq", "LtrUniqueIDstringPerCpStr"))
 
  OneGenePerCpCs <- foreach(id=identifier.v, .combine="cbind", .inorder=TRUE,
                            .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    
    id.header.ind <- grep(x=LtrUniqueIDstringPerCpStr, pattern=paste0("_", id, "_"))
    id.header <- LtrUniqueIDstringPerCpStr[id.header.ind]
    id.uniqueID <- LtrUniqueIDstringPerCpStr[id.header.ind+1L]
    rm("LtrUniqueIDstringPerCpStr", "id.header.ind"); gc()
    
    # Verify the order of cp/strength entries (should be smallest to greatest)
    # TO FOLLOW
    
    id.uniqueID.len <- length(id.uniqueID)
    
    # Convert unique IDs to name2 (HUGO gene symbols)
    # Output is list of df(uniqueID-name2) per cp/cs score
    uniqueID.name2.lst <- sapply(X=1:id.uniqueID.len, simplify=FALSE, 
                                 FUN=function(uniqueIDstr){
                                   uniqueIDstr <- id.uniqueID[uniqueIDstr]
                                   uniqueIDs <- as.numeric( strsplit(x=uniqueIDstr, split=";")[[1]] )
                                   df <- data.frame(uniqueID=uniqueIDs, name2=annotable[uniqueIDs], stringsAsFactors=FALSE)
                                   if( any(duplicated(df[["uniqueID"]])) | any(duplicated(df[["name2"]])) ){
                                     stop("Checkpoint 1.")
                                   }
                                   df
                                 })
    
    rm("id.uniqueID")
    
    fin.uniqueID.name2.lst <- as.list(rep(NA, times=id.uniqueID.len))
    # greatest cp/cs score taken by default
    fin.uniqueID.name2.lst[[id.uniqueID.len]] <- uniqueID.name2.lst[[id.uniqueID.len]]
    
    ref.gene.pool <- uniqueID.name2.lst[[id.uniqueID.len]][,"name2"]
    
    for( i in (id.uniqueID.len-1L):1 ){
      # Remove genes in current df that is present in the dfs with higher cp/cs
      drop.ind <- na.omit(
        match( ref.gene.pool, 
               table=uniqueID.name2.lst[[i]][["name2"]] ) 
      )
      fin.uniqueID.name2.lst[[i]] <- uniqueID.name2.lst[[i]][-(drop.ind),]
      ref.gene.pool <- unique(
        c(ref.gene.pool, 
          fin.uniqueID.name2.lst[[i]][,"name2"])
      )
    }
    
    # Check indeed 1:1 correspondence between gene and cs/cp
    
    temp <- do.call("rbind", fin.uniqueID.name2.lst)
    if( any(duplicated(temp[["uniqueID"]])) | any(duplicated(temp[["name2"]])) ){
      stop("Checkpoint 2.")
    }
    
    rm("temp", "uniqueID.name2.lst"); gc()
    
    #
    
    # Turn uniqueID and name2 list to strings
    fin.uniqueID.name2.lst <- lapply(X=fin.uniqueID.name2.lst, FUN=function(df){
      c(paste(df[[1]], collapse=";"), paste(df[[2]], collapse=";"))
    })
    
    # header - uniqueID - name2
    rbind( id.header, do.call("cbind", fin.uniqueID.name2.lst) )
    
  }
  
  write(x=unlist(OneGenePerCpCs[1:2,], use.names=FALSE), 
        file=paste0(output.dir, "/", gcb, "_", refseq, 
                    "_uniqueID_OneGenePerCpCs"))
  
  write(x=unlist(OneGenePerCpCs[c(1,3),], use.names=FALSE), 
        file=paste0(output.dir, "/", gcb, "_", refseq, 
                    "_name2_OneGenePerCpCs"))
}

# rm(list=ls())

  

 
  
  
  
  
  