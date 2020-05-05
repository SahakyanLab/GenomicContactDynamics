################################################################################
# KEGG enrichment analysis using enrichKEGG() from clusterProfiler package
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# BiocManager::install("clusterProfiler")
# library(clusterProfiler)

# clusterProfiler v3.8.1  For help: 
# https://yulab-smu.github.io/clusterProfiler-book/chapter1.html
# https://guangchuangyu.github.io/software/clusterProfiler
# Citation:
#  Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. 
#  clusterProfiler: an R package for comparing biological themes among gene clusters. 
#  OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
# KEGG supported organism listed in 
#  http://www.genome.jp/kegg/catalog/org_list.html’
### FUNCTION ###################################################################
enKEGG <- function(genes = genelist,
                   org = "hsa",
                   # c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")
                   inputKey = "ncbi-geneid",
                   # Background genes
                   bgr = uni,
                   # For saving table; if NULL, table not saved
                   filePath = "lib/ekegg.csv"){
  
  acceptedKeys <- c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")
  if(!inputKey%in%acceptedKeys){
    stop("Key type not accepted by KEGG. Choose among ",
         paste(acceptedKeys, collapse=","), ".")
  }
  
  ekegg <- enrichKEGG(gene=genes,
                      organism=org,
                      keyType=inputKey,
                      universe=bgr,
                      pAdjustMethod="BH",
                      pvalueCutoff=0.05, minGSSize=10, 
                      maxGSSize=500, qvalueCutoff=0.2,
                      #FALSE: download the latest KEGG data 
                      #For enrichment analysis
                      #TRUE: use KEGG.db
                      use_internal_data=FALSE)
  
  ekegg <- as.data.frame(ekegg)
  ekegg$Description <- gsub(x=ekegg$Description, pattern=",", replacement=" ", 
                            fixed=TRUE)
  
  # Save table
  if( !is.null(filePath[1]) ){
    write.csv(x=ekegg, file=filePath, row.names=FALSE, quote=FALSE)
  }
  
  return(ekegg)

}
