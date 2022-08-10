################################################################################
# Combine GO_ALL and KEGG enrichment analyses for a set of genes in one 
# dataframe.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(lib, "/funxAnno.R"))
# library(org.Hs.eg.db)
# library(clusterProfiler)
### FUNCTION ###################################################################
funxAnnoWrapper <- function(input = 'list(HUGO foreground genes, HUGO background genes)',
                            filePath = 'dir for saving individual GO and KEGG tables; 
                                        if NULL, tables not saved',
                            hugoEntrez.file = 'Conversion table of HUGO gene 
                                               symbols to ncbi-geneid for KEGG analysis'
){
  
  # Read conversion table for ncbi-geneid to HUGO symbols for output
  HE.conv <- read.delim(file=hugoEntrez.file, header=T, row.names=NULL, stringsAsFactors=F)
  HE.conv <- HE.conv[!is.na(HE.conv$SYMBOL) & !is.na(HE.conv$ENTREZID),]
  HE.conv$ENTREZID <- as.character(HE.conv$ENTREZID)
  
  # Gene sets
  GENE <- list()
  
  GENE$GO_ALL <- input
  GENE$GO_ALL <- lapply(X=GENE$GO_ALL, FUN=function(set){
    return( set[!is.na(set)] )
  }) 
  
  GENE$KEGG <- lapply(X=GENE$GO_ALL, FUN=function(set){
    set <- HE.conv$ENTREZID[ HE.conv$SYMBOL%in%set ]
    return( set[!is.na(set)] )
  })
  
  # Enrichment analysis
  funxAnnoOut <- sapply(X=names(GENE), simplify=F, FUN=function(app){
    
    chunk <- funxAnno(input=GENE[[app]],
                      org=ifelse(app=="GO_ALL", "org.Hs.eg.db", "hsa"), 
                      inputKey=ifelse(app=="GO_ALL", "SYMBOL", "ncbi-geneid"), 
                      approach=app, filePath=filePath)
    return(chunk)

  })
  
  funxAnnoOut <- do.call("rbind", funxAnnoOut)
  rownames(funxAnnoOut) <- NULL
  
  rm(GENE)
  
  # Convert ncbi geneids to HUGO
  KEGG.TF <- funxAnnoOut$ONTOLOGY=="KEGG"
  if( sum(KEGG.TF)>0 ){
    funxAnnoOut[KEGG.TF,"geneID"] <- sapply(X=funxAnnoOut[KEGG.TF, "geneID"], 
                                            simplify=TRUE, FUN=function(id){
                                              id <- unique(strsplit(x=id, split="/", fixed=TRUE)[[1]])
                                              paste(x=HE.conv[HE.conv$ENTREZID%in%id,"SYMBOL"],
                                                    collapse="/")
                                            })
  }
 
  return(funxAnnoOut)
  
} 
################################################################################

# rm(list=ls()); gc()