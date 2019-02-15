################################################################################ 
#Exploratory and other plots (SumStrengthOverNtis or meanvalsum)

################################################################################ 
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
#AnnotationVsPersistfiles 
ANNO.PERSIST.dir = paste0(objective.dir, "/out_AnnotationVsPersist")
#ANNO.PERSIST.dir = paste0(objective.dir, "/out_AnnotationVsPersist_sample")
output.dir = paste0(objective.dir, "/out_GeneVsNtis") 

library(foreach)
library(data.table)
library(dplyr)

#genome version
genome.ver = "hg19"

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify 
#a contact, two points should be far enough to filter for contacts within a TAD
gc.v = c("2", "05")

chr.v = "ALL" #c(1:22, "X", "ALL") 

data.v = c("ALL", "NM", "NR")

#check in AnnotationVsPersist files
#minimum overlap used when intersecting annotations to persistent contacts
olap.min = 1L 

#gc="2"
#chr="21"
#data="ALL"

################################################################################
################################################################################
#to convert gene names to entrezgene Ids
entrezId.convtbl <- fread(file=paste0(objective.dir, "/out_addEntrezidToAnno/",
                                      genome.ver, "anno_uniqueHGNC_entrezId.txt"),
                          data.table=FALSE, header=TRUE)
#persistence score
#specify this if you know the relevant ntis,
#if not get it from the AnnotationVsPersist files 
#should be characters
ntis.v <- as.character(1:21)

foreach(chr=chr.v) %do% {
  for(gc in gc.v){
    for(data in data.v){
      load(paste0(ANNO.PERSIST.dir,"/", "chr", chr, "_min", gc, "Mb_", 
                  genome.ver, "_", data, "_AnnoVsPersist_olapMin",
                  olap.min, ".RData"))
      ntis <- as.character(ANNO.PERSIST.MX[,"ntis"])
      
      #get unique ntis NOT present per annotation to save memory
      splt <- strsplit(x=ntis, split=",")
      presentntis <- lapply(X=splt, FUN=function(x){
        ntis.v[!ntis.v%in%as.numeric(x)]
      })
      
      num.rows <- nrow(ANNO.PERSIST.MX)
      if(length(presentntis)!=num.rows){
        stop("Checkpoint 1: Lengths don't match.")
      }
    
      #get genes per ntis 
      genesPerNtis <- list()
      for(ntis in ntis.v){
        ind <- sapply(X=presentntis, FUN=function(x)ntis%in%x)
        if(length(ind)==0){
          next
        } else {
          if(length(ind)!=num.rows){
            stop("Checkpoint 2: Lengths don't match.")
          }
        }
        genesPerNtis[[ntis]] <- unique( as.character(ANNO.PERSIST.MX[!ind,"name2"]) )
        #save genes per bin 
        write(genesPerNtis[[ntis]], 
              file=paste0(output.dir, "/", genome.ver, "_min", gc, 
                          "Mb_chr", chr, "_", data, "_PScore", ntis, 
                          "_GeneVsNtis.txt"))
      }
      #save corresponding files with entrezgeneIds
      for(ntis in names(genesPerNtis)){
        ind <- match(genesPerNtis[[ntis]], 
                     entrezId.convtbl[,"hgnc_symbol"])
        #do not remove genes with no entrezgene Id
        mx <- cbind( Gene=genesPerNtis[[ntis]], 
                     entrezgene=entrezId.convtbl[ind,"entrezgene"] )
        write.table(mx,
                    file=paste0(output.dir, "/entrezId_", genome.ver, 
                                "_min", gc, "Mb_chr", chr, "_", data, 
                                "_PScore", ntis, "_GeneVsNtis.txt"),
                    col.names=TRUE, row.names=FALSE, quote=FALSE)
        #code if using merge instead of match
        #gene.mx <- cbind(filler=rep("a"), Gene=genesPerNtis[[ntis]])
        #mx <- merge(gene.mx, entrezId.convtbl, 
        #            by.x="Gene", by.y="hgnc_symbol", 
        #            all.x=TRUE) #keep rows of x
        #mx <- set(mx, j="filler", value = NULL)
      }
  
    }
  }
}
  
#rm(list=ls())






