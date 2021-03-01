################################################################################
# Convert HUGO names to other gene identifiers 
# Set to Homo sapiens but can be applied to other organisms by changing the 
# dataset in useMart().
# Returns a table with origKeyType-desiredKeyType following the order and frquency
# of genes supplied
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#library(data.table)
#library(biomaRt)
### FUNCTION ###################################################################
convertGeneKeyType <- function(genes = genelist, 
                               origKeyType = "hgnc_symbol",  
                               # check biomaRt for proper reference to gene identifiers
                               convertTo = "entrezgene"){
  
  mrt <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                 dataset="hsapiens_gene_ensembl",
                 verbose=FALSE)
  
  output <- getBM(attributes = c(origKeyType, convertTo), 
                  filters = origKeyType, values = unique(genes), 
                  bmHeader = T, mart = mrt)
  # output: col1(filter/origKeyType), col2(convertTo)
  df <- by(output[,2], output[,1], 
           function(x) paste(x, collapse=","))
  df <- cbind(names(df), df)
  rownames(df) <- NULL
  colnames(df) <- c(origKeyType, convertTo)

  # make sure all genes in genelist are in table and same order
  genes <- data.frame(genes, stringsAsFactors=FALSE)
  colnames(genes) <- origKeyType
  df <- merge(x=genes, y=df, by=origKeyType, all.x=TRUE)
  df <- df[match(df[[origKeyType]], genes[[origKeyType]]),]
  # to change factors to characters
  return( rapply(object=df, f=as.character, 
               classes="factor", how="replace") )
}

# Other gene identifiers:
# uniprotswissprot
# ensembl_transcript_id

