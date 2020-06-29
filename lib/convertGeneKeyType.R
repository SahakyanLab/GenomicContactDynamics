################################################################################
# Convert HUGO names to other gene identifiers 
# Set to Homo sapiens but can be applied to other organisms by changing the 
# dataset in useMart().
# Returns a table with origKeyType-desiredKeyType following the order and frquency
# of genes supplied
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(org.Hs.eg.db)
# library(clusterProfiler)
### FUNCTION ###################################################################
convertGeneKeyType <- function(genes = NULL, # If genes=NULL & convTablePath specified, 
                               # whole conversion table returned.
                               origKeyType = "SYMBOL",  
                               convertTo = "ENTREZID",
                               org = "org.Hs.eg.db", # "hsa" (bitr_KEGG())
                               # If convTablePath=NULL, use bitr() or bitr_KEGG() 
                               # from clusterProfiler
                               convTablePath = NULL, #"lib/hg19anno_SYMBOLtoENTREZID_052020",
                               # If useKEGGAPI=TRUE, use bit_kegg()
                               useKEGGAPI = F,
                               drop.NA = F){
  
  if( !is.null(convTablePath) ){
    print("Using an existing conversion table...", quote=F)
    df <- read.delim(file=convTablePath, header=T, row.names=NULL, stringsAsFactors=F)
    if(drop.NA){
      df <- df[!is.na(df[,origKeyType]) & !is.na(df[,convertTo]),] 
    } 
    if( !is.null(genes) ){
      df <- df[df[,origKeyType]%in%genes,c(origKeyType, convertTo)]
      df <- df[order(match(x=df[,origKeyType], table=genes)),]
    } else {
      print("Full table loaded.", quote=F)
    }
    
  } else {
    print("Using bitr...", quote=F)
    
    if(useKEGGAPI==F){
      df <- bitr(geneID=genes, fromType=origKeyType, toType=convertTo, 
                 OrgDb=org, drop=drop.NA)
    } else {
      df <- bitr_kegg(geneID=genes, fromType=origKeyType, toType=convertTo, 
                      organism=org, drop=drop.NA)
    }
    if( any(duplicated(df[,origKeyType])) ){
      print("Multiple results for input...", quote=F)
      #df <- by(data=df[,convertTo], INDICES=df[,origKeyType], FUN=function(v){
      #  paste(v, collapse=",")
      #})
      #df <- stack(df)[,c("ind", "values")]
      #df$ind <- as.character(df$ind)
      #colnames(df) <- c(origKeyType, convertTo)
    }
  } 
  
  if(drop.NA){
    print("Dropping input if not converted...", quote=F)
    df <- df[!is.na(df[,1]) & !is.na(df[,2]),] 
  } 
  if( !identical(unique(df[,origKeyType]), genes) & drop.NA==F ){
    stop("Different order of genes in output.")
  } 
  
  rownames(df) <- NULL
  return(df)
  
}

#> keytypes(org.Hs.eg.db)
#[1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
#[9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
#[17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
#[25] "UNIGENE"      "UNIPROT"     

#> bitr_KEGG
#‘kegg’, ‘ncbi-geneid’ (same as ENTREZID), ‘ncbi-proteinid’ or ‘uniprot’
