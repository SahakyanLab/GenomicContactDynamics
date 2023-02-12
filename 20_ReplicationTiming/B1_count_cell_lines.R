################################################################################
# Count number of unique cell lines contributing RT data and if how many are tumors 
# and nontumors. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/20_ReplicationTiming")
norep.celldata.path = paste0(data.dir, "/replication_timing/out_process_datasets/hg19/norep.celldata.csv")
norep.tumor.path = paste0(data.dir, "/replication_timing/out_process_datasets/hg19/accessions.kept.norep.tumor.tsv")
norep.nontumor.path = paste0(data.dir, "/replication_timing/out_process_datasets/hg19/accessions.kept.norep.nontumor.tsv")
norep.allExcl.path = paste0(data.dir, "/replication_timing/out_process_datasets/hg19/accessions.kept.norep.allExcl.tsv")
norep.all.path = paste0(data.dir, "/replication_timing/out_process_datasets/hg19/accessions.kept.norep.all.tsv")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(stringr)
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
norep.celldata.df <- t(read.csv(norep.celldata.path, header=F, stringsAsFactors=F))
norep.tumor.df <- read.csv(norep.tumor.path, header=F, stringsAsFactors=F)
norep.nontumor.df <- read.csv(norep.nontumor.path, header=F, stringsAsFactors=F)
norep.allExcl.df <- read.csv(norep.allExcl.path, header=F, stringsAsFactors=F)
norep.all.df <- read.csv(norep.all.path, header=F, stringsAsFactors=F)

# Number of unique cell lines

cell.nmes.uniq <- sort(unique(tolower(
  str_replace_all(as.character(norep.celldata.df[,2]), " ", "_")
)))

# Number of unique tumor cell lines

# From A1_process_datasets.R
tumor.ct <- c("Patient Leukemia Sample", "Epithelial cells from colon carcinoma", 
              "Cervical Carcinoma", "Hepatocellular Carcinoma", "Leukemia",
              "Adenocarcinoma", "Neuroblastoma")

if( all(norep.tumor.df$V1 %in% norep.celldata.df[,1]) ){
  
  ind.celldata <- match(norep.tumor.df$V1, table=norep.celldata.df[,1])
  norep.tumor.df$cellline <- norep.celldata.df[ind.celldata,2]
  
  t.cell.nmes.uniq <- sort(unique(tolower(
    str_replace_all(as.character(norep.tumor.df[,"cellline"]), " ", "_")
  )))
  
  tumor.ct.uniq <- sort(unique(tolower(
    str_replace_all(as.character(tumor.ct), " ", "_")
  )))
  
  if( !identical(t.cell.nmes.uniq,tumor.ct.uniq) ){
    stop("Error with tumor cell line names.")
  }
  
  length(t.cell.nmes.uniq)
  
} 

# Number of unique non-tumor cell lines

if( all(norep.nontumor.df$V1 %in% norep.celldata.df[,1]) ){
  
  ind.celldata <- match(norep.nontumor.df$V1, table=norep.celldata.df[,1])
  norep.nontumor.df$cellline <- norep.celldata.df[ind.celldata,2]
  
  nt.cell.nmes.uniq <- sort(unique(tolower(
    str_replace_all(as.character(norep.nontumor.df[,"cellline"]), " ", "_")
  )))
  
  length(nt.cell.nmes.uniq)

} 

## These 5 cell lines were not classified into tumors or non tumor, see README
#exclude.ct <- c("LiverD16", "LiverD5", "LiverD8", "PancD12", "Epithelial cells immortalized with hTERT")



# rm(list=ls()); gc()