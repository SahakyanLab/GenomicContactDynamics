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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/16_GeneVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
genelist.dir = paste0(wk.dir, "/out_geneFilePerCp")
unmapped.dir = paste0(wk.dir, "/out_DAVID/unmapped")
out.dir = paste0(wk.dir, "/out_sampleGeneList")
### OTHER SETTINGS #############################################################
# gcb + refseq
prefix = "min2Mb_LTr_ALL"
SEED = 754 #587 #287
# Limit of genes that can be inputted to DAVID is 3000 but used 2999 cause 
#sometimes DAVID maps 3000 HUGO names to 3001 DAVID IDs.
N = 2999 
name.type = "name2" # "entrezID" | "name2"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
file.v <- list.files(path=genelist.dir, pattern=prefix, recursive=F)
file.v <- file.v[grep(x=file.v, pattern=name.type)]

for(fle in file.v){
  
  fle <- strsplit(x=fle, split=".", fixed=TRUE)[[1]][1]
  genes <- readLines(con=paste0(genelist.dir, "/", fle, ".txt"))

  # Save original gene list
  writeLines(unique(genes), con=paste0(out.dir, "/", fle, "_original.txt"))
  
  # If file of unmapped genes by DAVID is there, remove those genes.
  unmapped.file <- paste0(unmapped.dir, "/", fle, "_DAVIDunmapped.txt")
  if( file.exists(unmapped.file) ){
    
    print(paste0(fle, ": Removing unmapped then sampling."))
    
    unmapped <- readLines(con=unmapped.file)
    genes <- genes[!genes%in%unmapped]
    
    # Sample
    set.seed(SEED)
    sampled <- sample(x=genes, size=N, replace=FALSE)
    writeLines(sampled, con=paste0(out.dir, "/", fle, "_seed", SEED, "_",
                                   N, "sampled.txt"))
    
    rm(sampled)
    
  }
  
  print(fle, quote=FALSE)
  
  rm(unmapped.file, fle, genes)
  
}

# name 2; Cp=21 - out of 4221 --> 3941 DAVID IDs
# entrezID; Cp=21 - out of 3796 --> 3796 DAVID IDs

# rm(list=ls()); gc()