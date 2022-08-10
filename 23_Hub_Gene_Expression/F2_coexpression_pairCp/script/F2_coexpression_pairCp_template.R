################################################################################
# Get maximum Cp of all possible gene pairs brought together by a contact.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/CoreGenomeExplorer")

persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_coexpression_pairCp")

gene.id = "LTr_ALL" #"ALL"
annofilePath = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
bin.len = 40000
chr = "chrCHRREPLACE" #paste0("chr", c(1:22, "X"))
gcb = "min2Mb"
nCPU = 14 # Number of pairs per chr
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(gcb=="min2Mb"){
  min.gap <- 50
} else if(gcb=="min05Mb"){
  min.gap <- 12.5
} else {
  stop("No min.gap specified for gcb argument.")
}

out.name <- paste0(gcb, "_", gene.id, "_binlen", bin.len)

anno.df <- read.delim(file=annofilePath, header=T, stringsAsFactors=F, 
                      sep="\t")[,c("chrom", "name2", "txStart", "txEnd")]
if( any(duplicated(anno.df$name2)) ){
  stop("Duplicated genes in anno.df")
}

bin.start <- ceiling(anno.df$txStart/bin.len)
bin.end <- ceiling(anno.df$txEnd/bin.len)

#for(chr in chr.v){
  
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  
  ij.mx <- cbind(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
  rm(PERSIST.MX)
  gc()
  
  testi <- ij.mx[,"i"]%in%c(bin.start, bin.end)
  testj <- ij.mx[,"j"]%in%c(bin.start, bin.end)
  
  ij.mx <- ij.mx[testi & testj,]
  rm(testi, testj)
  
  # All indices wrt to anno.df
  chr.g.ind.v <- which(anno.df$chrom==chr)
  genes.len <- length(chr.g.ind.v)
  rm(anno.df)
  
  pair.ind.df <- expand.grid(chr.g.ind.v, chr.g.ind.v, stringsAsFactors=F)
  pair.ind.df <- pair.ind.df[pair.ind.df$Var1 < pair.ind.df$Var2,]
  pair.len <- length(pair.ind.df$Var1)
  
  if( pair.len != (genes.len*genes.len-genes.len)/2 ){
    stop("Number of gene pairs not equal to expected.")
  }
  
  print(paste0(chr, ": ", pair.len, " pairs."), quote=F)
  
  rm(chr.g.ind.v)
  gc()
  
  toExport <- c("pair.ind.df", "ij.mx", "bin.start", "bin.end")
  
  #### PARALLEL EXECUTION #########
  
  PAIRCP.MX <- foreach(pair.ind.v=isplitVector(1:pair.len, chunks=nCPU), 
                       .inorder=F, .combine="rbind",
                       .export=toExport, .noexport=ls()[!ls()%in%toExport]
                        
  ) %op% {
    
    chunk <- sapply(X=pair.ind.v, simplify=F, FUN=function(pair.ind){
      
      print(pair.ind, quote=F)
      
      g1.ind <- pair.ind.df[pair.ind,1]
      g2.ind <- pair.ind.df[pair.ind,2]
      
      g1.bins <- bin.start[g1.ind]:bin.end[g1.ind]
      g2.bins <- bin.start[g2.ind]:bin.end[g2.ind]
      
      max.gap.genes <- diff( range(c(g1.bins, g2.bins)) )
      
      if(max.gap.genes <= min.gap){
        
        maxcp <- NA
        
      } else {
        
        test.g12 <- ij.mx[,"i"]%in%g1.bins & ij.mx[,"j"]%in%g2.bins
        test.g21 <- ij.mx[,"i"]%in%g2.bins & ij.mx[,"j"]%in%g1.bins
        testij <- test.g12 | test.g21
        rm(test.g12, test.g21)
        
        maxcp <- ifelse( sum(testij)==0,  NA, 
                         max(as.numeric(ij.mx[testij,"Cp"]), na.rm=F) )
        
      }
      
      return( c(g1.ind, g2.ind, maxcp) )
      
    })
    
    return(do.call("rbind", chunk))
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  print("Done with PAIRCOR.MX parallel execution.", quote=F)
  
  dimnames(PAIRCP.MX)[[1]] <- NULL
  dimnames(PAIRCP.MX)[[2]] <- c("gene1.ind.anno", "gene2.ind.anno", "maxCp")
  
  rm(pair.ind.df, ij.mx, bin.start, bin.end)
  gc()
  
  PAIRCP.MX <- PAIRCP.MX[order(PAIRCP.MX[,1],PAIRCP.MX[,2]), 1:length(PAIRCP.MX[1,])]
  
  if( any(PAIRCP.MX[,2] <= PAIRCP.MX[,1]) | pair.len != length(PAIRCP.MX[,1]) ){
    stop("Error in final object.")
  } else {
    save(PAIRCP.MX, file=paste0(out.dir, "/", chr, "_", out.name, "_pairMaxCp_coexpression.RData"))
  }

#} # chr.v for loop end

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()