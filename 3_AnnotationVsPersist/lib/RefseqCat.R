#A RefSeq record has a distinct accession number format that begins with two
#characters followed by an underscore (e.g., NP_)
#These characters indicate the RefSeq category of the record/annotation.
#For more info, check the RefSeq FAQ of NCBI

#The script reports the RefSeq categories (accession prefixes) present in the
#annotation file

################################################################################
#install.packages("data.table") 
library(data.table) #for fread and rbindlist

################################################################################
################################################################################
RefSeqCat <- function(filePath = "/Users/ltamon/Database/ucsc_tables/hg19anno",
                      genomeVersion = "hg19",
         #accessionCol refers to name of column containing the accession numbers
                      accessionCol = "name" ){ 
  anno.file <- fread(file=filePath, header=TRUE, data.table=FALSE)
  
  anno.type <- unique( sapply( as.vector(anno.file[,accessionCol]), 
                               simplify=TRUE, USE.NAMES=FALSE, 
                               function(x){
                                 b <- strsplit(x, split="_")
                                 strsplit(unlist(b), split="\t")[[1]]
                               } ) )
  
  anno.type <- paste(anno.type, collapse=", ")
  print(paste0("RefSeq categories in ", genomeVersion, " annotation file are ", 
               anno.type, "."))
}
