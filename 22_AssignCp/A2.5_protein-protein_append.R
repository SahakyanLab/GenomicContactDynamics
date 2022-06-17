################################################################################
# Append maximum Cp data to FunCoup full data format
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Avoid left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/X1_")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/X1_")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
data.dir = paste0(home.dir, "/Database")
full.file = paste0(data.dir, "/FunCoup/FC5.0_H.sapiens_full")
append.file = paste0(wk.dir, "/out_protein-protein/min2Mb_protein-protein_maxCp.csv")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for( x in c("full", "append") ){
  
  eval(parse(text=paste0( "file.nme <- ", x, ".file")))
  df <- data.table::fread(file=file.nme, stringsAsFactors=F, header=T, data.table=F, check.names=T)
  
  if( !all(grepl(c(df$X2.Gene1, df$X3.Gene2), pattern="ENSG")) ){
    rm(df); stop(paste0(x, ".df: Checkpoint 1."))
  } 
  
  df$X2.Gene1 <- as.numeric(gsub(df$X2.Gene1, pattern="ENSG", replacement=""))
  df$X3.Gene2 <- as.numeric(gsub(df$X3.Gene2, pattern="ENSG", replacement=""))
  id <- paste0(df$X2.Gene1, "-", df$X3.Gene2)
  if( any(duplicated(id)) ){
    rm(df); stop(paste0(x, ".df: Checkpoint 2."))
  } 
  eval(parse(text=paste0( x, ".id <- id")))
  eval(parse(text=paste0( x, ".df <- df")))
  
  rm(df, id, file.nme)
  
}

colnmes <- c("X4.Gene1", "X5.Gene2", "chrom", "Cp")
filler <- matrix(data=NA, nrow=length(full.df[,1]), ncol=length(colnmes), dimnames=list(NULL, colnmes))
full.df <- cbind.data.frame(filler, full.df)
rm(filler)
full.df[match(x=append.id, table=full.id), 1:4] <- append.df[,colnmes]

out.file <- gsub(x=append.file, pattern=".csv", replacement="_full.csv") 
write.csv(full.df, file=out.file, row.names=F)

# rm(list=ls()); gc()