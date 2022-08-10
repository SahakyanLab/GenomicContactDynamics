################################################################################
# Pairwise correlation of cross-tissue expression values of all possible gene 
# pairs in expression data. Use correlation as a measure of coexpression. 
# percNDopenUpperLimPair argument is an open upper limit for allowable fraction
# of tissues with no expression value for the gene pair. If percNDopenUpperLimPair
# equal to 0.2, only gene pairs with no expression values < 0.2 of tissues 
# are selected. Expression values < expr.cutoff converted to 0. 
# Output is a unique genes x unique genes matrix of correlation values with 
# only the lower triangle filled. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

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

exprData.dir = paste0(wk.dir, "/out_cleanExprData")
out.dir = paste0(wk.dir, "/out_coexpression_pairCor")

gene.id = "LTr_ALL" #"ALL"
annofilePath = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
src.id = "data2" 
expr.cutoff = 0
nCPU = 1 #4 # Number of pairs # 4G per core
cor.meth = "pearson" # "pearson" "spearman"
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
                   cor.meth)

exprDataPath <- paste0(exprData.dir, "/expr_", src.id ,"_cutoff0_", gene.id, ".csv")
expr.df <- read.csv(file=exprDataPath, header=T, stringsAsFactors=F)
expr.df$chr <- NULL

anno.df <- read.delim(file=annofilePath, header=T, stringsAsFactors=F, sep="\t")[,"name2", drop=F]
#anno.df <- anno.df[20928:20929, "name2", drop=F]
expr.df <- expr.df[expr.df$Gene.Name%in%anno.df$name2,]

# Append to expr.df missing genes from anno.df

toaddgenes <- setdiff(anno.df$name2, expr.df$Gene.Name)
appendgenes.df <- as.data.frame(matrix(data=NA, nrow=length(toaddgenes), 
                                       ncol=ncol(expr.df),
                                       dimnames=list(NULL, colnames(expr.df)))
)
appendgenes.df$Gene.Name <- toaddgenes

expr.df <- rbind.data.frame(expr.df, appendgenes.df, stringsAsFactors=F)
if( all(anno.df$name2%in%expr.df$Gene.Name) ){
  expr.df <- expr.df[match(x=anno.df$name2, table=expr.df$Gene.Name),]
}

if( !identical(anno.df$name2, expr.df$Gene.Name) ){
  stop("Not identical order of genes in expr.df and anno.df.")
}

rm(appendgenes.df, toaddgenes, anno.df, exprDataPath)

# Convert to matrix

tiss.TF <- !colnames(expr.df)%in%c("chr", "Gene.Name")
expr.mx <- data.matrix(expr.df[,tiss.TF])
dimnames(expr.mx)[[1]] <- expr.df$Gene.Name

# 
expr.mx[expr.mx < expr.cutoff] <- 0

rm(expr.df, tiss.TF)
gc()

tiss.len <- length(expr.mx[1,])
expr.genes.len <- length(expr.mx[,1])

toExport <- c("expr.mx", "cor.meth", "expr.genes.len", "tiss.len")

#### PARALLEL EXECUTION #########

PAIRCOR.MX <- foreach(g1.ind.v=isplitVector(1:(expr.genes.len-1), chunks=nCPU), 
                      .inorder=F, .combine="rbind",
                      .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  chunk <- sapply(X=g1.ind.v, simplify=F, FUN=function(g1.ind){
    
    print(g1.ind, quote=F)
    
    g2.ind.v <- (g1.ind+1):expr.genes.len
  
    chunk1 <- foreach(g2.ind=g2.ind.v, .combine="rbind"
    
    ) %do% {
      
      frTisspairAtleast1NA <- colSums(expr.mx[c(g1.ind, g2.ind),])
      frTisspairAtleast1NA <- sum(is.na(frTisspairAtleast1NA), na.rm=F)/tiss.len
  
      cor.val <- cor(x=expr.mx[g1.ind,], y=expr.mx[g2.ind,], method=cor.meth)
      
      return( c( g1.ind, g2.ind, cor.val, frTisspairAtleast1NA) )
    
    }
    
    return(chunk1)
    
  })
  
  return(do.call("rbind", chunk))
  
}

### END OF PARALLEL EXECUTION ###

print("Done with PAIRCOR.MX parallel execution.", quote=F)

# Store in matrix (upper triangle only)
dimnames(PAIRCOR.MX)[[1]] <- NULL
dimnames(PAIRCOR.MX)[[2]] <- c("g1.ind", "g2.ind", "cor", "frTisspairAtleast1NA")
PAIRCOR.MX <- PAIRCOR.MX[order(PAIRCOR.MX[,1],PAIRCOR.MX[,2]), 1:length(PAIRCOR.MX[1,]), 
                         drop=F]

if( any(PAIRCOR.MX[,2] <= PAIRCOR.MX[,1]) ){
  stop("Invalid indices.")
}

#PAIRCOR.MX <- matrix(data=NA, nrow=expr.genes.len, ncol=expr.genes.len, byrow=F,
#                     dimnames=list( dimnames(expr.mx)[[1]], 
#                                    dimnames(expr.mx)[[1]] )
#                     )
#PAIRCOR.MX[lower.tri(PAIRCOR.MX)] <- lowertri.val

NAcorpair.len <- sum(is.na(PAIRCOR.MX[,"cor"]))
pair.len <- length(PAIRCOR.MX[,1])
stat.df <- data.frame(NAcorpair=NAcorpair.len, totpairInAnno=pair.len, 
                      frNAcorpair=NAcorpair.len/pair.len, stringsAsFactors=F)
write.table(x=stat.df, file=paste0(out.dir, "/", out.name, "_NApair.txt"), 
            row.names=F, quote=F)

rm(expr.mx, NAcorpair.len, pair.len, stat.df)
gc()

save(PAIRCOR.MX, file=paste0(out.dir, "/", out.name, "_coexpression.RData"))

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()
