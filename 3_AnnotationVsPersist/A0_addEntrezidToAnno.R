################################################################################
#Add entrezgene ID (NCBI Entrez Gene database) to the annotation files (from UCSC)

################################################################################
#functions to source
lib <- "/Users/ltamon/ProjectBoard/lib"
#annotation files
annoFile.path <- "/Users/ltamon/Database/ucsc_tables"
output.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist/out_addEntrezidToAnno"

library(data.table)

################################################################################
source(paste0(lib, "/convertHGNCsymb.R"))

################################################################################
################################################################################
for(anno.nme in c("hg19anno", "hg38anno")){
  anno.file <- fread(file=paste0(annoFile.path, "/", anno.nme), 
                     header=TRUE, data.table=FALSE)
  genelist <- sort(unique(anno.file[,"name2"])) 
  
  #save unique HGNC symbols from annotation file
  write(genelist, file=paste0(output.dir, "/", anno.nme, 
                              "_uniqueHGNC.txt"))

  #make table of unique HGNC and corresponding entrezgeneid
  hgnc.entrezid <- convertHGNCsymb(genes = genelist, 
                                   convertTo = c("entrezgene"))
  
  hgnc.entrezid.uniq <- by(hgnc.entrezid[[1]][,"entrezgene"], 
                           hgnc.entrezid[[1]][,"hgnc_symbol"], 
                           function(x) paste(x, collapse = ","))
  tbl <- cbind(names(hgnc.entrezid.uniq), hgnc.entrezid.uniq)
  colnames(tbl) <- c("hgnc_symbol", "entrezgene")
  write.table(tbl, file=paste0(output.dir, "/", anno.nme, 
                               "_uniqueHGNC_entrezID.txt"),
              sep="\t", quote=FALSE,
              col.names=TRUE, row.names=FALSE)
  
  #add entrezgeneid to annotation file 
  anno.entrezid <- merge(anno.file, as.data.frame(hgnc.entrezid[[1]]), 
                         by.x="name2", by.y="hgnc_symbol", 
                         all.x=TRUE) #keep rows of original annotation file
  write.table(anno.entrezid, file=paste0(output.dir, "/", anno.nme, 
                                         "_entrezID.txt"),
              sep="\t", quote=FALSE,
              col.names=TRUE, row.names=FALSE)
}

#rm(list=ls())

