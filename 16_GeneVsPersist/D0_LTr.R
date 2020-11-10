################################################################################
# Identify longest transcript per gene
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    annofile.dir = "/t1-data/user/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
### OTHER SETTINGS #############################################################
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
refseq.v = c("ALL", "NM", "NR")
anno.nme = "hg19anno"
### FUNCTION ###################################################################
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(refseq in refseq.v){
  annotable <- fread(file=paste0(annofile.dir, "/", anno.nme, "_", refseq), 
                     header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)
  annotable$txSize <- annotable$txEnd-annotable$txStart
 
  LtrPerHugo <- by(data=annotable[, c("uniqueID", "txSize", "name")], 
                   INDICES=annotable$name2, 
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
  
  annotable <- annotable[annotable$uniqueID%in%LtrPerHugo,]
  
  write.table(x=annotable, file=paste0(annofile.dir, "/", anno.nme, "LTr_", 
                                       refseq),
              sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

# rm(list=ls())

  



