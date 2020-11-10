################################################################################
# Examine UCSC annotation table
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables"
    os = "Mac"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}

output.dir = annofile.dir
### OTHER SETTINGS #############################################################
annofile.v = "hg19anno" 
anno.nme = "hg19anno"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(anno.nme in annofile.v){
  
  # HUGO symbols both NM and NR
  # hg19: 37559
  HUGO.NM <- fread(file=paste0(annofile.dir, "/", anno.nme, "_NM"), 
                     header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)[["name2"]] 
  # hg19: 9335
  HUGO.NR <- fread(file=paste0(annofile.dir, "/", anno.nme, "_NR"), 
                   header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)[["name2"]] 
  
  # hg19: 1478; intersect returns unique
  HUGO.olap <- intersect(unique(HUGO.NM), unique(HUGO.NR))
  cat(length(HUGO.olap), "HUGO symbols overlapping.")
  
}


anno.file <- fread(file=paste0(annofile.dir, "/", anno.nme, "_ALL"), 
                 header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)

# Annotation table column names
c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
  "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", 
  "cdsEndStat", "exonFrames", "uniqueID")
