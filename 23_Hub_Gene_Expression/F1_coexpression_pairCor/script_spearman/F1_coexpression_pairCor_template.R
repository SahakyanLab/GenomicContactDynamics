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

exprData.dir = paste0(wk.dir, "/out_cleanExprData")
out.dir = paste0(wk.dir, "/out_coexpression_pairCor")

gene.id = "LTr_ALL" #"ALL"
annofilePath = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
chr = "chrCHRREPLACE"
cor.meth = "spearman" # "pearson" "spearman"
expr.cutoff = 0 #0.5
nCPU = 1 #4 # Number of pairs # 4G per core
src.id = "data2" 
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
                   chr, "_", cor.meth)

exprDataPath <- paste0(exprData.dir, "/expr_", src.id ,"_cutoff0_", gene.id, ".csv")
expr.df <- read.csv(file=exprDataPath, header=T, stringsAsFactors=F)
expr.df$chr <- NULL

anno.df <- read.delim(file=annofilePath, header=T, stringsAsFactors=F, sep="\t")[,c("chrom", "name2")]
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
} else {
  expr.chr <- anno.df$chrom
}

rm(appendgenes.df, toaddgenes, anno.df, exprDataPath)

# Convert expr.df to matrix

tiss.TF <- !colnames(expr.df)%in%c("chr", "Gene.Name")
expr.mx <- data.matrix(expr.df[,tiss.TF])
dimnames(expr.mx)[[1]] <- expr.df$Gene.Name

# 
expr.mx[expr.mx < expr.cutoff] <- 0

#
tiss.len <- sum(tiss.TF)

rm(expr.df, tiss.TF)
gc()

# Generate gene pairs per chr

chr.g.ind.v <- which(expr.chr==chr)
genes.len <- length(chr.g.ind.v)

pair.ind.df <- expand.grid(chr.g.ind.v, chr.g.ind.v, stringsAsFactors=F)
pair.ind.df <- pair.ind.df[pair.ind.df$Var1 < pair.ind.df$Var2,]
pair.len <- length(pair.ind.df$Var1)

if( pair.len != (genes.len*genes.len-genes.len)/2 ){
  stop("Number of gene pairs not equal to expected.")
}

print(paste0(chr, ": ", pair.len, " pairs."), quote=F)

rm(chr.g.ind.v, expr.chr)
gc()

toExport <- c("pair.ind.df", "expr.mx", "tiss.len", "cor.meth")

#### PARALLEL EXECUTION #########

PAIRCOR.MX <- foreach(pair.ind.v=isplitVector(1:pair.len, chunks=nCPU), 
                      .inorder=F, .combine="rbind",
                      .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  chunk <- sapply(X=pair.ind.v, simplify=F, FUN=function(pair.ind){
    
    print(pair.ind, quote=F)
    
    g1.ind <- pair.ind.df[pair.ind,1]
    g2.ind <- pair.ind.df[pair.ind,2]
    
    frTisspairAtleast1NA <- colSums(expr.mx[c(g1.ind, g2.ind),], na.rm=F)
    noNAtisspair.TF <- !is.na(frTisspairAtleast1NA)
    frTisspairAtleast1NA <- sum(!noNAtisspair.TF, na.rm=F)/tiss.len
    
    cor.val <- cor(x=expr.mx[g1.ind, noNAtisspair.TF], 
                   y=expr.mx[g2.ind, noNAtisspair.TF], 
                   method=cor.meth)
    
    return( c( g1.ind, g2.ind, cor.val, frTisspairAtleast1NA) )
  
  })
  
  return(do.call("rbind", chunk))
  
}

### END OF PARALLEL EXECUTION ###

print("Done with PAIRCOR.MX parallel execution.", quote=F)

dimnames(PAIRCOR.MX)[[1]] <- NULL
dimnames(PAIRCOR.MX)[[2]] <- c("g1.ind", "g2.ind", "cor", "frTisspairAtleast1NA")
PAIRCOR.MX <- PAIRCOR.MX[order(PAIRCOR.MX[,1],PAIRCOR.MX[,2]), 1:length(PAIRCOR.MX[1,])]

rm(pair.ind.df)
gc()

# 

NAcorpair.len <- sum(is.na(PAIRCOR.MX[,"cor"]))
stat.df <- data.frame(NAcorpair=NAcorpair.len, totpairInAnno=pair.len, 
                      frNAcorpair=NAcorpair.len/pair.len, stringsAsFactors=F)
write.table(x=stat.df, file=paste0(out.dir, "/", out.name, "_NApair.txt"), 
            row.names=F, quote=F)

rm(expr.mx, stat.df)
gc()

# cor can be NA if either x or y have sd = 0
#if( !identical(is.na(PAIRCOR.MX[,"cor"]),
#                     PAIRCOR.MX[,"frTisspairAtleast1NA"]==1) ){
# stop("Error in cor or frTisspairAtleast1NA calculation.")
#}

if( pair.len != length(PAIRCOR.MX[,1]) |
    any(PAIRCOR.MX[,"g2.ind"] <= PAIRCOR.MX[,"g1.ind"]) |
    any(is.na(PAIRCOR.MX[,"frTisspairAtleast1NA"])) 
  ){
  stop("Error in final object.")
} else {
  save(PAIRCOR.MX, file=paste0(out.dir, "/", out.name, "_coexpression.RData"))
}

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()



