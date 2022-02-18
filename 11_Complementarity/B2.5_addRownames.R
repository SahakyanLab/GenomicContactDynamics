################################################################################
# Add rownames to CII.MX align
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
compl.dir = out.dir = paste0(wk.dir, "/out_constraints")
# File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = "chr17" #paste0("chr", c(1:22, "X"))
type = "align"
bin.len = 40000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Chromosome length file
chrLen.df <- read.table(file=chrLenfile, as.is=TRUE, header=TRUE)

for(chr in chr.v){
  
  # CII.MX kmer as reference
  load(file=paste0(compl.dir, "/", chr, "_kmer_", gcb, ".RData"))
  rwK <- rownames(CII.MX)
  rm(CII.MX);gc()
  
  # Load CII.MX
  load(file=paste0(compl.dir, "/", chr, "_", type, "_", gcb, ".RData"))
  
  tot.bin <- ceiling(chrLen.df[chrLen.df$chromosome==chr, "length.bp"]/bin.len)
  ubins <- 1:tot.bin
  
  # Create matrix of all pairs 
  contact.mx <- data.matrix(expand.grid(ubins, ubins))
  dimnames(contact.mx)[[1]] <- 1:nrow(contact.mx)
  # To remove self pairs and duplicates like {1,2}(keep) and {2,1}
  contact.mx <- contact.mx[contact.mx[,2]-contact.mx[,1]>0,]
  
  if( nrow(contact.mx)==nrow(CII.MX) ){
    
    dimnames(CII.MX)[[1]] <- rownames(contact.mx)
    rm(tot.bin, ubins, contact.mx); gc()
    
    if(!identical(rownames(CII.MX), rwK)){
      stop("Checkpoint 1.")
    } else {
      save(CII.MX, file=paste0(out.dir, "/", chr, "_", type, "_", gcb, ".RData"))
    }
    
  } else {
    
    stop("Checkpoint 1.")
    
  }
  
  rm(CII.MX); gc()
  print(paste0(chr, " done!"), quote=FALSE)
  
}

# rm(list=ls()); gc()
