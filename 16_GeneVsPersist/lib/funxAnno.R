################################################################################
# Functional analyses of genes via GO or KEGG
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(lib, "/enKEGG.R"))
# source(paste0(lib, "/enGO.R"))
### FUNCTION ###################################################################
funxAnno <- function(input = "list of foregound and background genes",
                     # c("SYMBOL", "kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")
                     inputKey = "SYMBOL",
                     approach = "GO_BP",
                     # For saving GO or KEGG table; if NULL, table not saved
                     filePath = "lib/ekegg.txt"
){
  
  # Determine and apply requested functional analysis method 
  if( grepl(pattern="KEGG", x=approach) ){
    
    if(inputKey!="ncbi-geneid"){ "Approach KEGG but input key not ncbi-geneid." }
    # KEGG enrichment analysis   
    funxAnnoOut <- enKEGG(genes=input[[1]],
                          org="hsa",
                          # refers to the type of input or the keyType in enrichKEGG
                          # one of "kegg", "ncbi-geneid", "ncib-proteinid" and "uniprot"
                          inputKey="ncbi-geneid",
                          # for saving table
                          filePath=filePath,
                          bgr=input[[2]])
                                      
  } else if( grepl(pattern="GO", x=approach) ){
    
    # GO enrichment analysis
    funxAnnoOut <- enGO(genes=input[[1]],
                        # Organism Db object
                        orgDb=org.Hs.eg.db,
                        # Refers to the type of input or the keyType in enrichGO
                        # Is it Entrez ID? HUGO symbols? etc.
                        inputKey="SYMBOL",
                        # Vector of GO Domains/subontologies of interest
                        ontology=gsub(pattern="GO_", replacement="", 
                                      x=approach, fixed=TRUE),
                        # Background genes
                        bgr=input[[2]],
                        # For saving table
                        filePath=filePath)
                            
  } else {
    stop("Invalid approach supplied. Choose between GO or KEGG.")
  }
  
  return(as.data.frame(funxAnnoOut, stringsAsFactors=FALSE))
  
}
################################################################################
