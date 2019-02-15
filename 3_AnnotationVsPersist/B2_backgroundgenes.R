################################################################################
#make list of background genes (with HiC21 data) for functional annotation

################################################################################
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist/"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
ANNO.PERSIST.dir = paste0(objective.dir, "/out_AnnotationVsPersist")
#ANNO.PERSIST.dir = paste0(objective.dir, "/out_AnnotationVsPersist_sample")
#output.dir = paste0(objective.dir, "/out_backgroundgenes_test")
output.dir = paste0(objective.dir, "/out_backgroundgenes")

library(data.table) #for fread
library(foreach)

#genome version of annotation file used to derive gene list
genome.ver = "hg19"

#2(2MB gap) or "05"(0.5 MB gap), refers to minimum gap between contacting bins, 
#two points should be far enough to filter for contacts within a TAD
gc.v = c("2", "05")

chr.v = "ALL" 

data.v = c("ALL", "NM", "NR")

#check AnnotationVsPersist filenames
olap.min = 1L 

################################################################################
#for converting to entrezgene IDs
#to convert gene names to entrezgene Ids
entrezId.convtbl <- fread(file=paste0(objective.dir, "/out_addEntrezidToAnno/",
                                      genome.ver, "anno_uniqueHGNC_entrezId.txt"),
                          data.table=FALSE, header=TRUE)
for(gc in gc.v){
  for(chr in chr.v){
    foreach(data=data.v, .inorder=TRUE) %do% {
      load(file=paste0(ANNO.PERSIST.dir,"/", "chr", chr, "_min", gc, "Mb_", 
                       genome.ver, "_", data, "_AnnoVsPersist_olapMin", 
                       olap.min, ".RData"))
        genes <- unique(ANNO.PERSIST.MX[,"name2"])
        write(genes, file=paste0(output.dir, "/", genome.ver, "_min", gc, 
                         "Mb_chr", chr, "_", data, "_background.txt"))
        #store background genes with entrezId
        ind <- match(genes, entrezId.convtbl[,"hgnc_symbol"])
        #do not remove genes with no entrezgene Id
        mx <- cbind(Gene=genes, 
                    entrezgene=entrezId.convtbl[ind,"entrezgene"])
        write.table(mx, file=paste0(output.dir, "/entrezId_", genome.ver, "_min", gc, 
                                    "Mb_chr", chr, "_", data, "_background.txt"),
                    col.names=TRUE, row.names=FALSE, quote=FALSE)
        #code if using merge instead of match
        #gene.mx <- cbind(filler=rep("a"), Gene=genes)
        #mx <- merge(gene.mx, entrezId.convtbl, 
        #            by.x="Gene", by.y="hgnc_symbol", 
        #            all.x=TRUE) #keep rows of x
        #mx <- set(mx, j="filler", value = NULL)
    }
  }
}

#rm(list=ls())

  

