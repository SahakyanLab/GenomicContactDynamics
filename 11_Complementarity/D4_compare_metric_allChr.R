################################################################################
# Correlate C|| metrics combining contacts from all chromosomes.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/11_Complementarity"
  } else if(whorunsit == "LiezelCluster"){
    #lib = "/t1-data/user/ltamon/DPhil/lib"
    #data.dir = "/t1-data/user/ltamon/Database"
    #wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/pending/11_Constraints"
    lib = "/stopgap/sahakyanlab/ltamon/DPhil/lib"
    data.dir = "/stopgap/sahakyanlab/ltamon/Database"
    wk.dir = "/stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/pending/11_Constraints"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
compl.dir = paste0(wk.dir, "/z_ignore_git/out_constraints/merged_final")
#compl.dir = paste0(wk.dir, "/out_constraints/merged_final")
out.dir = paste0(wk.dir, "/out_compare_metric_allChr")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr0.v = unique(paste("chr", c(1:22, "X"), sep="")) # Shuffled later on
nCPU = 1 # chr and comb.v (n=9) below 
plotOnly = TRUE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(plotOnly==FALSE){
  
  set.seed(1234)
  chr.v <- sample(x=chr0.v, size=length(chr0.v), replace=FALSE)
  chr.v.len <- length(chr.v)
  
  if( any(duplicated(chr.v)) | any(!chr0.v%in%chr.v) ){
    stop("Checkpoint 1.")
  }
  
  toExport <- c("chr.v", "compl.dir", "gcb")
  
  DF <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                .inorder=FALSE, .combine="rbind",
                .export=toExport, .noexport=ls()[!ls()%in%toExport]
                
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      chr <- chr.v[i]
      
      # Load CII.MX kmer
      load(file=paste0(compl.dir, "/", chr, "_kmer_", gcb, ".RData"))
      
      # Keep Cp column to know which contacts are short-lange or long-range
      df <- CII.MX[, c("Cp", "C||", "Gfree")]
      colnames(df) <- c("Cp", "kmer", "Gfree")
      rm(CII.MX)
      gc()
      
      # Total contacts, same for CII.MX align
      len0 <- length(df[,1])
      
      # Load CII.MX align
      load(file=paste0(compl.dir, "/", chr, "_align_", gcb, ".RData"))
      if( len0!=length(CII.MX[,1]) ){
        stop(paste0(chr, ": kmer and align CII.MX have different dimensions."))
      }
      if( !identical(as.numeric(df[,"Cp"]), as.numeric(CII.MX[,"Cp"])) ){
        stop(paste0(chr, ": kmer and align CII.MX have different Cp values."))
      }
      df <- cbind(df, align=CII.MX[,"C||"])
      rownames(df) <- NULL
      rm(CII.MX, len0)
      gc()
      
      # Take only contacts with values for 3 complementarity metrics
      df <- df[!is.na(df[,"kmer"]) & !is.na(df[,"Gfree"]) & !is.na(df[,"align"]),]
      
      return(df)
      
    })
    
    return( do.call("rbind", chunk) )
    
  }
  
  save(DF, file=paste0(out.dir, "/", gcb, "_CIImetricValues.RData"))
  
  LR.TF <- !is.na(DF[,"Cp"])
  SR.TF <- is.na(DF[,"Cp"])
  nonNA.len <- length(DF[,1])
  DF <- as.data.frame(DF[,colnames(DF)!="Cp"])
  
  # Calculate correlations
  
  comb.v <- combn(x=c("kmer", "align", "Gfree"), m=2, FUN=paste, collapse="-")
  comb.v.len <- length(comb.v)
  ijset <- c("All", "SR", "LR")
  
  comb.v <- paste(rep(x=ijset, each=comb.v.len), rep(x=comb.v, times=length(ijset)), sep="-")
  comb.v.len <- length(comb.v)
  
  toExport <- c("comb.v", "DF", "nonNA.len", "LR.TF", "SR.TF")
  
  COR <- foreach(itr=isplitVector(1:comb.v.len, chunks=nCPU), 
                 .inorder=FALSE, .combine="rbind",
                 .export=toExport, .noexport=ls()[!ls()%in%toExport]
                 
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      comb.id <- comb.v[i]
      comb <- strsplit(x=comb.id, split="-", fixed=TRUE)[[1]]
      
      if(comb[[1]]=="All"){
        ij.TF <- rep(x=TRUE, times=nonNA.len)
      } else {
        eval(parse(text=paste0( 'ij.TF <- ', comb[[1]], '.TF' )))
      }
      
      eval(parse(text=paste0("lm.obj <- lm(", comb[[2]], "~", comb[[3]], ", data=DF[ij.TF,])")))
      coef.v <- setNames(object=coef(lm.obj)[c("(Intercept)", comb[[3]])], nm=c("intercept", "slope"))
      
      return(
        
        data.frame(`ijset.y.x`=comb.id, intercept=coef.v[["intercept"]], slope=coef.v[["slope"]], 
                   Rsquared=summary(lm.obj)$r.squared,
                   pearson=cor(x=DF[ ij.TF,comb[[3]] ], y=DF[ ij.TF,comb[[2]] ], method="pearson"),
                   #kendall=cor(x=DF[ ij.TF,comb[[3]] ], y=DF[ ij.TF,comb[[2]] ], method="kendall"),
                   spearman=cor(x=DF[ ij.TF,comb[[3]] ], y=DF[ ij.TF,comb[[2]] ], method="spearman"),
                   stringsAsFactors=FALSE) 
        
      )
      
    })
    
    return( do.call("rbind", chunk) )
    
  }
  
  save(COR, file=paste0(out.dir, "/", gcb, "_CIImetricPairwiseCor.RData"))
  
  
} else {
  
  load(file=paste0(out.dir, "/", gcb, "_CIImetricPairwiseCor.RData"))
  
}

# Plot
df <- reshape2::melt(COR, measure.vars=c("pearson", "spearman"))
tmp <- do.call("rbind", strsplit(x=df$ijset.y.x, split="-", fixed=TRUE))
colnames(tmp) <- c("ijset", "y", "x")

df <- cbind(tmp, df, yx=paste(tmp[,"y"], tmp[,"x"], sep="-"))
df$ijset <- factor(x=df$ijset, levels=c("All", "LR", "SR"))

#ggplot(data=df, aes(x=yx, y=abs(value))) +
ggplot(data=df, aes(x=yx, y=value)) +
  geom_col(aes(fill=ijset), position="dodge", colour="white") + 
  scale_fill_manual(values=c("#F39B7FFF", "#00A087FF", "#8491B4FF")) + 
  #labs(x=NULL, fill=NULL, title=paste0(gcb, "_AbsoluteOfcorrelationsTaken_valueforGfreeNeg_chr1To22andX")) +
  labs(x=NULL, fill=NULL, title=paste0(gcb, "_chr1To22andX")) +
  bgr2 +
  theme(strip.background=element_blank(),
        axis.text.x=element_text(angle=0, size=15)) + 
  facet_grid(.~variable) 

ggsave(filename=paste0(out.dir, "/", gcb, "_CIImetricPairwiseCor.pdf"),
       width=10, height=10, unit="in")
# rm(list=ls()); gc()
