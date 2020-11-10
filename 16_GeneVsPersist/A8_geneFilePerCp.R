################################################################################
# Get list of genes colocalising with Cp of interest
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
genelist.dir = paste0(wk.dir, "/out_anno_union")
out.dir = paste0(wk.dir, "/out_geneFilePerCp")
### OTHER SETTINGS #############################################################
# gcb + refseq
prefix = "min2Mb_ALL"
# Forebackground strings defining foreground and background, foreground;background
# to differentiate TLC from LC 
cp.str.v = c("cp_19", "cp_20", "cp_21") #"cp_21" #c("cp_19", "cp_20", "cp_21") 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cp.str.v <- unique(cp.str.v)
cp.str.v.len <- length(cp.str.v)
pat.v <- paste0( "_", unlist(strsplit(x=cp.str.v, split=";")), "_end")
flenme <- paste(cp.str.v, collapse="_")
for(type in c("name2", "entrezID", "uniqueID")){
  genes.str <- readLines(con=paste0(genelist.dir, "/", prefix, "_", type))
  ind <- sapply(X=pat.v, FUN=grep, x=genes.str, fixed=TRUE)
  if( length(ind)!=cp.str.v.len ){ stop("Checkpoint 1.") }
  
  genes <- list()
  for(i in ind){
    genes[[as.character(i)]] <- strsplit(x=genes.str[i+1L], split=";")[[1]] 
  }
  genes <- sort( unique(unlist(genes, use.names=FALSE)) )
  writeLines(genes, con=paste0(out.dir, "/", prefix, "_", flenme, "_", type, ".txt"))
  
  rm(genes.str, genes, ind); gc()
  print(type, quote=FALSE)
}

# rm(list=ls())