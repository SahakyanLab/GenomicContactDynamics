################################################################################
# Generate a matrix of coexpression (measured via pearson/spearman correlation) 
# of all possible gene pairs in expression data.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=TRUE)

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

exprData.dir = paste0(wk.dir, "/out_cleanExprData")
out.dir = paste0(wk.dir, "/out_coexpression_pairwise")

gene.id = "LTr_ALL" #"ALL"
annofilePath = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
src.id = "data2" 
expr.cutoff = 0
nCPU = 1 # Number of pairs
cor.meth = "pearson" # "pearson" "spearman"
percNDopenUpperLimPair = 1
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
out.name <- paste0("expr_", src.id, "_cutoff", expr.cutoff, "_", gene.id, "_",
                   cor.meth, "_percNDopenUpperLimPair", percNDopenUpperLimPair)

exprDataPath <- paste0(exprData.dir, "/expr_", src.id ,"_cutoff0_", gene.id, ".csv")
expr.df <- read.csv(file=exprDataPath, header=T, stringsAsFactors=F)[1:10,]

anno.df <- read.delim(file=annofilePath, header=T, stringsAsFactors=F, sep="\t")[,c("chrom", "name2")]
expr.df <- expr.df[expr.df$Gene.Name%in%anno.df$name2,]
expr.df <- expr.df[order(expr.df$Gene.Name, decreasing=F),]

tiss.TF <- !colnames(expr.df)%in%c("chr", "Gene.Name")
expr.mx <- data.matrix(expr.df[,tiss.TF])
dimnames(expr.mx)[[1]] <- expr.df$Gene.Name

expr.mx[expr.mx < expr.cutoff] <- 0

rm(anno.df, expr.df, tiss.TF)
gc()

tiss.len <- length(expr.mx[1,])
expr.genes.len <- length(expr.mx[,1])

toExport <- c("expr.mx", "cor.meth", "expr.genes.len", "tiss.len")

#### PARALLEL EXECUTION #########

PAIRCOR.MX <- foreach(p1.v=isplitVector(1:(expr.genes.len-1), chunks=nCPU), 
                      .inorder=F, .combine="rbind",
                      .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  chunk <- sapply(X=p1.v, simplify=F, FUN=function(p1){
    
    pair.ind <- (p1+1):expr.genes.len
  
    chunk1 <- foreach(p2=pair.ind, .combine="rbind"
    
    ) %do% {
      
      NAtisspair <- colSums(expr.mx[c(p1,p2),])
      NAtisspair <- sum(is.na(NAtisspair), na.rm=F)/tiss.len
  
      if( NAtisspair < percNDopenUpperLimPair ){
        cor.val <- cor(x=expr.mx[p1,], y=expr.mx[p2,], method=cor.meth)
      } else {
        cor.val <- NA
      }
      
      return( c( p1, p2, cor.val) )
    
    }
    
    return(chunk1)
    
  })
  
  return(do.call("rbind", chunk))
  
}

### END OF PARALLEL EXECUTION ###

# Store in matrix (upper triangle only)
dimnames(PAIRCOR.MX)[[1]] <- NULL
lowertri.val <- PAIRCOR.MX[order(PAIRCOR.MX[,1],PAIRCOR.MX[,2]),3]
PAIRCOR.MX <- matrix(data=NA, nrow=expr.genes.len, ncol=expr.genes.len, byrow=F,
                     dimnames=list( dimnames(expr.mx)[[1]], 
                                    dimnames(expr.mx)[[1]] )
                     )
PAIRCOR.MX[lower.tri(PAIRCOR.MX)] <- lowertri.val

NApair.len <- sum(is.na(lowertri.val))
pair.len <- length(lowertri.val)
stat.df <- data.frame(NApair=NApair.len, totpairInExpr=pair.len, 
                      frNApair=NApair.len/pair.len, stringsAsFactors=F)
write.table(x=stat.df, file=paste0(out.dir, "/", out.name, "_NApair.txt"), 
            row.names=F, quote=F)

rm(expr.mx, lowertri.val, NApair.len, pair.len, stat.df)
gc()

save(PAIRCOR.MX, file=paste0(out.dir, "/", out.name, "_coexpression.RData"))

# rm(list=ls()); gc()
