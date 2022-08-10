################################################################################
# Identify whether gene pairs are within a persistent hub or not
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
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

# Use PAIRCOR.MX to get gene pair indices (pearson, data2 version)
pairCor.dir = paste0(wk.dir, "/out_coexpression_pairCor") 
out.dir = paste0(wk.dir, "/out_coexpression_pairHub1")

gene.id = "LTr_ALL" 
hubsum.dir = paste0(wk.dir, "/out_hubsummary/All_LTr")
annofilePath = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
bin.len = 40000
chr.v = "chr20" #paste0("chr", c(1:22, "X"))
gcb = "min2Mb"
hub.id = "min2Mb_All_topCP3_gapBin50"
nCPU = 1 # Number of pairs per chr
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
out.name <- paste0("HUBID", hub.id, "_", gene.id)

gene.v <- read.delim(file=annofilePath, header=T, stringsAsFactors=F, sep="\t")[,c("name2")]
if( any(duplicated(gene.v)) ){
  stop("Duplicated genes.")
}

hsum.df <- read.csv(file=paste0(hubsum.dir, "/", hub.id, "_hubsum.csv"), 
                    stringsAsFactors=F, header=T)

for(chr in chr.v){
  
  load(paste0(pairCor.dir, "/expr_data2_cutoff0.5_", gene.id, "_", chr, 
              "_pearson_coexpression.RData"))
  
  pair.ind.df <- PAIRCOR.MX[,c("g1.ind", "g2.ind")]
  pair.len <- length(pair.ind.df[,1])
  
  gPerHub <- hsum.df[hsum.df$chr==chr,"genes"]
  
  if( length(gPerHub)==0 ){
    
    print(paste0(chr, ": No hubs, skipped."), quote=F)
    next
    
  } else {
    
    gPerHub <- strsplit(x=gPerHub, split=";", fixed=T)
    
    rm(PAIRCOR.MX)
    gc()
    
    toExport <- c("pair.ind.df", "gPerHub", "gene.v")
    
    #### PARALLEL EXECUTION #########
    
    PAIRHUB.MX <- foreach(pair.ind.v=isplitVector(1:pair.len, chunks=nCPU), 
                          .inorder=F, .combine="rbind",
                          .export=toExport, .noexport=ls()[!ls()%in%toExport]
                          
    ) %op% {
      
      chunk <- sapply(X=pair.ind.v, simplify=F, FUN=function(pair.ind){
        
        #print(pair.ind, quote=F)
        
        g1.ind <- pair.ind.df[pair.ind,1]
        g2.ind <- pair.ind.df[pair.ind,2]
        
        g1 <- gene.v[g1.ind]
        g2 <- gene.v[g2.ind]
        
        if(g1 == g2){
          stop(paste0(chr, " Pair ", pair.ind, ": Same gene1 and gene2."))
          rm(gPerHub)
        }
        
        hub.TF <- lapply(X=gPerHub, FUN=function(hub){
          c(g1%in%hub, g2%in%hub)
        })
        hub.TF <- do.call("rbind", hub.TF)
        hub.TF <- any(rowSums(hub.TF)==2)
        
        return( c(g1.ind, g2.ind, as.numeric(hub.TF)) )
        
      })
      
      return(do.call("rbind", chunk))
      
    }
    
    ### END OF PARALLEL EXECUTION ###
    
    dimnames(PAIRHUB.MX) <- list(NULL, c("g1.ind", "g2.ind", "WithinHub.TF"))
    
    rm(pair.ind.df, gPerHub)
    gc()
    
    PAIRHUB.MX <- PAIRHUB.MX[order(PAIRHUB.MX[,1],PAIRHUB.MX[,2]), 1:length(PAIRHUB.MX[1,])]
    
    if( any(PAIRHUB.MX[,2] <= PAIRHUB.MX[,1]) | any(is.na(PAIRHUB.MX[,"WithinHub.TF"])) ){
      stop("Error in PAIRHUB.MX")
    } else {
      save(PAIRHUB.MX, file=paste0(out.dir, "/", chr, "_", out.name, "_pairHub_coexpression.RData"))
    }
    
    rm(PAIRHUB.MX)
    gc()
    
    print(paste0(chr, " done!"), quote=F)
    
  }

} # chr.v for loop end

# rm(list=ls()); gc()
