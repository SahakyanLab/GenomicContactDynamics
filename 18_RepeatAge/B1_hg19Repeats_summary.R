################################################################################
# Convert repeat masker file to .RData for subsequent analyses.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/3_RepeatAge"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
repmaskfile = paste0(data.dir, "/ucsc_tables/hsa_RepeatMasker/RepeatMasker_hg19")
out.dir = paste0(wk.dir, "/out_hg19Repeats_summary")
### OTHER SETTINGS #############################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
repmasker.df <- fread(file=repmaskfile, header=TRUE, data.table=FALSE, 
                      stringsAsFactors=FALSE)

# to separate the same repNames with different repFamily and repClass 
# repClass repFamily repName
# Simple_repeat	Simple_repeat	(CATTC)n
# Satellite	Satellite	(CATTC)n
# Simple_repeat	Simple_repeat	(GAATG)n
# Satellite	Satellite	(GAATG)n
repmasker.df$uniqueID <- apply(X=repmasker.df[,c("repClass", "repFamily", "repName")],
                               MARGIN=1, FUN=function(x)paste(x, collapse=";"))

# reformat table - unique repName(subfamily)
# use uniqueID as indices
col.v <- setdiff(x=colnames(repmasker.df), 
                 y="uniqueID")
lst <- sapply(X=col.v, simplify=FALSE, FUN=function(col){
  out <- by(data=repmasker.df[,col], 
            INDICES=repmasker.df[,"uniqueID"], 
            FUN=function(x){
              
              if(col=="bin"){
                c(bin=paste(x, collapse=";"),
                  copyNumber=length(x))
              } else if(col%in%c("repClass", "repFamily", "repName")){
                # there are duplicated repName with different class and family
                paste(unique(x), collapse=";")
              } else {
                paste(x, collapse=";")
              }
              
            })
  
  if(col=="bin"){
    out <- do.call(rbind, out)
  } else {
    as.vector(out)
  }
  
})

# nrow=1397 but 1395 unique repNames
REPEAT.MX <- do.call("cbind.data.frame", lst)
setnames(x=REPEAT.MX, old="bin.bin", new="bin")
setnames(x=REPEAT.MX, old="bin.copyNumber", new="copyNumber")
# unfactorize
REPEAT.MX <- rapply(REPEAT.MX, as.character, classes="factor", how="replace")
save(REPEAT.MX, file=paste0(out.dir, "/hg19repeats_repName_base.RData"))

# rm(list=ls()); gc()

