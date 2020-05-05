################################################################################
# GO enrichment analysis using enrichGO() from clusterProfiler package
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(org.Hs.eg.db)
# BiocManager::install("clusterProfiler"); library(clusterProfiler)
# clusterProfiler v3.8.1  For help: 
# https://yulab-smu.github.io/clusterProfiler-book/chapter1.html
# https://guangchuangyu.github.io/software/clusterProfiler
# Citation:
#  Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. 
#  clusterProfiler: an R package for comparing biological themes among gene clusters. 
#  OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
### FUNCTION ###################################################################
enGO <- function(genes = genelist,
                 # Organism Db object
                 orgDb = org.Hs.eg.db,
                 inputKey = "SYMBOL", # "SYMBOL" | "ncbi-geneid"
                 #Domains/subontologies:
                 # Biological Process (BP)
                 # Cellular Component (CC)
                 # Molecular Function (MF)
                 ontology = "BP", # "BP" | "CC" | "MF"
                 # Background genes
                 bgr = uni,
                 # For saving table; if NULL, table not saved
                 filePath = "lib/ego.csv"){
  
  acceptedKeys <- c("SYMBOL", "ncbi-geneid")
  if(!inputKey%in%acceptedKeys){
    stop("Key type not accepted by KEGG. Choose among ",
         paste(acceptedKeys, collapse=","), ".")
  }
  inputKey <- ifelse(inputKey=="ncbi-geneid", "ENTREZID", inputKey)
  
  ego <- enrichGO(gene=unique(genes),
                  OrgDb=orgDb,
                  keyType=inputKey,
                  ont=ont,
                  pvalueCutoff=0.05,
                  pAdjustMethod="BH",
                  # Gene set size limits for the GO terms to be considered (DEF).
                  # This was added so that parent terms like "biological process"
                  # won't be tested.
                  minGSSize=10, 
                  maxGSSize=500,
                  universe=unique(bgr))
  
  ego <- data.frame(funxAnnoOut, stringsAsFactors=FALSE, row.names=NULL)
  ego$Description <- gsub(x=ego$Description, pattern=",", replacement=" ",
                          fixed=TRUE)
  # Save table
  if( !is.null(filePath) ){
    write.csv(x=ego, file=filePath, row.names=FALSE, quote=FALSE)
  } 
  return(ego)
    
}
