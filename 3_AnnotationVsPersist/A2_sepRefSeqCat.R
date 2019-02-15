################################################################################
#Separate AnnotationVsPersist file per RefSeq categories present

################################################################################ 
#functions to source
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
#annotation file 
#annoFile = "/Users/ltamon/Database/ucsc_tables/hg19anno"
annoFile = "/t1-home/icbml/ltamon/Database/ucsc_tables/hg19anno"
#AnnotationVsPersistfiles 
ANNO.PERSIST.dir = paste0(objective.dir, "/out_AnnotationVsPersist")
#ANNO.PERSIST.dir = paste0(objective.dir, "/out_AnnotationVsPersist_test")

library(foreach)
library(data.table) #for fread and rbindlist

#genome version
genome.ver <- "hg19"

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify 
#a contact, two points should be far enough to filter for contacts within a TAD
gc.v <- c("2", "05")

chr.v <- c("ALL") #c(1:22, "X", "ALL") 

#check AnnotationVsPersist filenames
olap.min <- 1L 

################################################################################
source(paste0(lib, "/RefSeqCat.R")) 
source(paste0(lib, "/UCSCTableReadFilter.R"))

################################################################################ 
#to check for RefSeq categories in the annotation file, 
#use the RefSeqCat function
RefSeqCat(filePath = annoFile,
          genomeVersion = genome.ver,
          accessionCol = "name")

#separate ANNO.PERSIST.MX based on RefSeq Categories
for(gc in gc.v){
  for(chr in chr.v){
    for(data in c("NM", "NR")){
      load(file=paste0(ANNO.PERSIST.dir,"/", "chr", chr, "_min", gc, "Mb_", 
                       genome.ver, "_ALL_AnnoVsPersist_olapMin", olap.min, ".RData"))
      if(data=="NM"){
        scheme <- "Human.Nuclear.mRNA"
      } else {
        scheme <- "Human.ncRNA"
      }
      ANNO.PERSIST.MX <- UCSCTableReadFilter(Table=ANNO.PERSIST.MX,
                                             Filtering.Scheme=scheme)
      save(ANNO.PERSIST.MX, file=paste0(ANNO.PERSIST.dir, "/chr", chr, "_min", gc, "Mb_", 
                                        genome.ver, "_", data, "_AnnoVsPersist_olapMin", 
                                        olap.min, ".RData"))
    }
  }
}

#rm(list=ls())




