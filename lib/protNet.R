#Potein network analysis using STRINGdb package

#using STRINGdb (database of known and predicted protein-protein interactions)
#The interactions include direct (physical) and indirect (functional) associations

################################################################################
#package for Installing and Managing Bioconductor Packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("STRINGdb", version = "3.8")
library(STRINGdb)

################################################################################
################################################################################
#VECTOR of genes as input
get.protnet <- function(genes = genelist,
                        outputfile = "lib/protNet.pdf"){
  
  df <- data.frame(genes=genes)
  #map the genes to the STRING database identifiers
  #$map(my_data_frame, my_data_frame_id_col_names, takeFirst = , 
  #removeUnmappedRows = , quiet = )
  #quiet=FALSE, issue warning about unmapped values
  str.mapped <- string_db$map( df, "genes", removeUnmappedRows=TRUE, quiet=FALSE )
  pdf(height=8, width=8, file=outputfile)
  #extract genes mapped (specify how many) and produce an image of the 
  #STRING network for those
  string_db$plot_network(
    str.mapped$STRING_id[1:min(400, length(str.mapped$STRING_id))]
  )
  dev.off()
}

#rm(list=ls())




