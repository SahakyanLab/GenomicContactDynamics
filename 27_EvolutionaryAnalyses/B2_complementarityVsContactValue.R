################################################################################
# Intersect complementarity values (CII.MX) with any data of format:
# chr-start.i-end.i-start.j-end.j-value
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/27_EvolutionaryAnalyses")
CII.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity/out_constraints_hg38_GfreeSingleNorm/merged_final")
value.dir = paste0(data.dir, "/Phylo-HMRF/out_combine_splitPerChr") # Path to a directory
out.dir = paste0(wk.dir, "/out_complementarityVsContactValue") 
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh38_chr_info.txt")
### OTHER SETTINGS #############################################################
value.id = "genome_state_Phylo-HMRF_contact50K_norm_KR"
coord.system = "zero-based"
gcb = "min2Mb"
compl.data.ids = "kmer"
compl.types = c("C||", "Gfree", "sdDifference")
chr.v = paste0("chr", c(17,19,21:22))
bin.len = 50000
value.header = T
options(scipen=999) # To prevent CII.MX and VALUE.df rownames to become sci note
# needed for comparing the two
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrlen.df <- read.table(file=chrlen.file, stringsAsFactors=F, header=T)

for(chr in chr.v){
  
  value.df <- fread(file=paste0(value.dir, "/", chr, "_", value.id, ".txt"), 
                    header=value.header, stringsAsFactors=F, data.table=F)
  
  #
  
  if( any(value.df[,1] != chr) ){
    rm(value.df)
    stop(paste0(chr, ": Mixture of chr."))
  } 
  
  #
  
  bin.len.check <- unique(c(unique(value.df[,3] - value.df[,2]), unique(value.df[,5] - value.df[,4])))
  
  if( !identical(as.numeric(bin.len.check), as.numeric(bin.len)) ){
    stop(paste0(chr, ": Wrong bin length."))
  } else if( coord.system == "zero-based" & bin.len.check == bin.len ){
    
    value.df[,2] <- value.df[,2] + 1
    value.df[,4] <- value.df[,4] + 1
  
  } else if( coord.system == "one-based" & (bin.len.check + 1) != bin.len ){
  } else {
    stop(paste0(chr, ": Wrong bin length or coord.system argument."))
  }
   
  # Convert coordinates to bin numbers
  
  value.df[,2:5] <- ceiling(value.df[,2:5] / bin.len)
  check.val <- unique(unique(c(value.df[,3] - value.df[,2])), unique(c(value.df[,5] - value.df[,4])))
  
  if(check.val != 0){
    stop(paste0(chr, ": Region map to multiple bin numbers."))
  } else{
    value.df[,1] <- value.df[,2] <- value.df[,4] <- NULL
  }
  
  #
  
  value.names <- setdiff(colnames(value.df), colnames(value.df)[1:2])
  colnames(value.df) <- c( "i", "j", value.names)
  
  #
  
  value.df <- value.df[value.df$i != value.df$j, ]
  print(paste0(chr, ": Self-pairs i == j removed."), quote=F)
  
  if( any(value.df$j < value.df$i) ){
    rm(value.df)
    stop(paste0(chr, ": i > j entries."))
  }
  
  # Merge value.df with VALUE.df keeping order of latter, which should match with CII.MX
  
  tot.bin <- ceiling( chrlen.df$length.bp[chrlen.df$chromosome == chr] / bin.len )
  mx <- matrix(data=NA, nrow=tot.bin, ncol=tot.bin, dimnames=list(i=1:tot.bin, j=1:tot.bin))
  
  rw.len <- length(value.df$i)
  value.df$value.index <- 1:length(value.df$i)
  for(rw in 1:rw.len){
    mx[ value.df$i[[rw]], value.df$j[[rw]] ] <- value.df$value.index[[rw]]
  }
  
  VALUE.df <- reshape2::melt(mx)
  rm(mx)
  VALUE.df <- VALUE.df[VALUE.df$j - VALUE.df$i > 0,] # To remove self pairs and duplicates like {1,2}(keep) and {2,1}
  VALUE.df$i <- as.numeric(VALUE.df$i)
  VALUE.df$j <- as.numeric(VALUE.df$j)
  
  # Actual merging
  
  VALUE.df$row.num <- as.numeric(rownames(VALUE.df))
  merged.tmp <- merge(x=VALUE.df, y=value.df[,setdiff(colnames(value.df), c("i", "j"))],
                      by.x="value", by.y="value.index", all.x=T)
  VALUE.df <- merged.tmp[ match(VALUE.df$row.num, table=merged.tmp$row.num), ]
  rownames(VALUE.df) <- as.character(VALUE.df$row.num)
  
  rm(merged.tmp, value.df)
  gc()
  
  #
  
  CII.tmp <- list()
  for(id in compl.data.ids){
    
    load(paste0(CII.dir, "/", chr, "_", id, "_", gcb, ".RData"))
    if( ! identical( as.matrix(VALUE.df[,c("i", "j")]), CII.MX[,c("i", "j")] ) ){
      rm(VALUE.df)
      stop(paste0(chr, " ", id, ": VALUE.df and CII.MX contact order not matching."))
    }
    dimnames(CII.MX)[[2]] <- paste0(id, ".", dimnames(CII.MX)[[2]])
    
    CII.tmp[[id]] <- CII.MX[ ,intersect(paste0(id, ".", compl.types), dimnames(CII.MX)[[2]]), drop=F ] 
    rm(CII.MX)
    
    print(paste0(chr, " ", id, ": Adding data."), quote=F)
    
  }
  
  names(CII.tmp) <- NULL
  CII.tmp <- do.call("cbind.data.frame", CII.tmp)
  VAL.MX <- cbind.data.frame(VALUE.df[,c("i", "j")], CII.tmp, 
                             VALUE.df[,value.names, drop=F])
  
  rm(CII.tmp, VALUE.df)
  gc()
  
  save(VAL.MX, file=paste0(out.dir, "/", chr, "_", gcb, "_", value.id, ".RData"))
  
  rm(VAL.MX)
  gc()
  
  print(paste0(chr, " done!"), quote=F)
  
}

# rm(list=ls()); gc()
