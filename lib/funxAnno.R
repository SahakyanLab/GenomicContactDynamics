################################################################################
# Functional analyses of genes via GO or KEGG (enrichGO and enrichKEGG is more 
# reliable with more updated data than many other tools. - Guangchuang Yu,
# clusterProfiler developer)
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
funxAnno <- function(input = list(foreground, background),
                     org = org.Hs.eg.db, # org.Hs.eg.db (GO, human) | "hsa" (KEGG, human)
                     # c("SYMBOL", "kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")
                     inputKey = "SYMBOL",
                     approach = "GO_BP",
                     # For saving GO or KEGG table; if NULL, table not saved
                     filePath = "lib/ekegg.txt"
){
  
  approach <- gsub(x=approach, pattern="GO_", replacement="", fixed=TRUE)
  accepted.KEGG <- c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot")
  accepted.BP <- accepted.CC <- accepted.MF <- accepted.ALL <- c("SYMBOL", "ncbi-geneid") 
  eval(parse(text=paste0(
    "acceptedKeys <- accepted.", approach
  )))
 
  if(!inputKey%in%acceptedKeys){
    stop("Key type not accepted by", approach, ". Choose among ",
         paste(acceptedKeys, collapse=", "), ".")
  }
  
  # Clean gene lists
  input <- lapply(X=input, FUN=function(genes){
    unique(genes[!is.na(genes) & genes!="NA" & genes!="" & genes!=" "])
  })
  
  if(approach=="KEGG"){
    
    e.out <- enrichKEGG(gene=input[[1]], organism=org, keyType=inputKey,
                        # clusterProfiler uses the same adjustment methods 
                        # as the R function stats::p.adjust 
                        # The "BH" (aka "fdr") and "BY" method of Benjamini, Hochberg, 
                        # and Yekutieli control the false discovery rate, the expected 
                        # proportion of false discoveries amongst the rejected hypotheses. 
                        # The false discovery rate is a less stringent condition than 
                        # the family-wise error rate, so these methods are more powerful 
                        # than the others.
                        pvalueCutoff=0.05, pAdjustMethod="BH", universe=input[[2]],
                        # Gene set size limits for the GO terms to be considered (DEF).
                        # This was added so that parent terms like "biological process"
                        # won't be tested. 
                        minGSSize=10, maxGSSize=500, qvalueCutoff=0.2,
                        #If FALSE, download the latest KEGG data; If TRUE, use KEGG.db
                        use_internal_data=FALSE)
   
  } else if( approach%in%c("BP", "CC", "MF", "ALL") ){
    
    e.out <- enrichGO(gene=input[[1]], OrgDb=org, 
                      keyType=ifelse(inputKey=="ncbi-geneid", "ENTREZID", inputKey), 
                      ont=approach, pvalueCutoff=0.05, pAdjustMethod="BH", 
                      universe=input[[2]], minGSSize=10, maxGSSize=500,
                      # Whether mapping gene ID to gene Name
                      readable=ifelse(inputKey=="ncbi-geneid", TRUE, FALSE),
                      # If ont='ALL', whether pool 3 GO sub-ontologies
                      pool=ifelse(approach=="ALL", TRUE, FALSE))
    
  } else {
    stop("Invalid approach supplied. Choices are GO_BP, GO_CC, GO_MF, GO_ALL and KEGG.")
  }
  
  e.out <- data.frame(e.out, stringsAsFactors=FALSE, row.names=NULL)
  # To avoid separation of description names when opening the csv file.
  e.out$Description <- gsub(x=e.out$Description, pattern=",", replacement=" ",  
                            fixed=TRUE)
  if(approach!="ALL" & nrow(e.out)>0){
    e.out <- cbind.data.frame(ONTOLOGY=approach, e.out, stringsAsFactors=FALSE)
  }
 
  # Save table
  if( !is.null(filePath) ){
    write.csv(x=e.out, file=filePath, row.names=FALSE, quote=FALSE)
  }
  
  return(e.out)
  
}
################################################################################
