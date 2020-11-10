################################################################################
# Sample 3000 genes from gene list for DAVID functional clustering
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
genelist.dir = paste0(wk.dir, "/out_geneFilePerCp")
unmapped.dir = paste0(wk.dir, "/out_DAVID")
out.dir = paste0(wk.dir, "/out_sampleGeneList")
### OTHER SETTINGS #############################################################
# gcb + refseq
prefix = "min2Mb_ALL"
SEED = 921
N = 3000
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
file.v <- list.files(path=genelist.dir, pattern=prefix)
file.v <- file.v[grep(x=file.v, pattern="entrezID")]

for(fle in file.v){
  fle <- strsplit(x=fle, split=".", fixed=TRUE)[[1]][1]
  genes <- readLines(con=paste0(genelist.dir, "/", fle, ".txt"))
  # Remove unmapped genes by DAVID
  unmapped <- read.delim(file=paste0(unmapped.dir, "/", fle, "_genesNotInOut.txt"),
                         stringsAsFactors=FALSE, header=TRUE, sep="\t")
  unmapped <- rownames(unmapped)
  genes <- genes[!genes%in%unmapped]
  set.seed(SEED)
  sampled <- sample(x=genes, size=N, replace=FALSE)
  writeLines(sampled, con=paste0(out.dir, "/", fle, "_seed", SEED, "_sampled.txt"))
  print(fle, quote=FALSE)
}

# rm(list=ls())