################################################################################
# Plot gene expression data vs. Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = paste0(wk.dir, "/out_geneExprVsCp")
exprData.dir = paste0(wk.dir, "/out_cleanExprData")
# Unique HUGO genes overlapping with Cp regions (any type of overlap)
#CpGenesPath = paste0(wk.dir, "/txt_CpGenes/min2Mb_ALL_name2")
CpGenes.dir = paste0(wk.dir, "/txt_CpGenes")
### OTHER SETTINGS #############################################################
gene.id = "LTr_ALL" #"ALL"
src.id = "data2" 
gcb = "min2Mb"
expr.cutoff = 0.5 # TPM/FPKM below this means not expressed
col.v = c(MEAN="#d9c741", MEDIAN="#edb340", SD="#32a864", SDdivMEAN="#32a8a4", 
          ND="gray50", NE="gray70", LE="#2171B5", ME="#9ECAE1", HE="#A50F15")
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
exprDataPath <- paste0(exprData.dir, "/expr_", src.id ,"_cutoff0_", gene.id, ".csv")
out.id <- paste0(gcb, "_", src.id, "_", gene.id, "_geneExprVsCp")

if(plotOnly==FALSE){
  
  temp <- readLines(con=paste0(CpGenes.dir, "/", gcb, "_", gene.id, "_name2"))
  CpGenesInd <- grep(x=temp, pattern=">all_genes_cp_", fixed=TRUE)
  # String of genes per Cp
  CpGenes <- temp[CpGenesInd + 1L]
  names(CpGenes) <- gsub(x=temp[CpGenesInd], pattern=">all_genes_cp_|_end", replacement="")
  rm(temp, CpGenesInd)
  gc()
  CpGenes.len <- length(CpGenes)
  # List of genes overlapping with all long-range contact regions and regions
  # of Cp 1 to 21. 
  CpGenes <- strsplit(x=CpGenes, split=";", fixed=TRUE)
  
  expr.df <- read.csv(file=exprDataPath, header=TRUE, stringsAsFactors=FALSE)
  tiss.TF <- !colnames(expr.df)%in%c("chr", "Gene.Name")
  
  # Classify genes in expr.df into Cp categories and calculate values
  CPGNEXPR.DF <- sapply(X=names(CpGenes), simplify=FALSE, FUN=function(cp){
    
    if( any(duplicated(CpGenes[[cp]])) ){
      stop(paste0("Duplicated genes in ", cp))
    }
    
    # No. of Cp genes
    n0 <- length(CpGenes[[cp]])
    # No. of Cp genes not in expr.df
    nND <- sum(!CpGenes[[cp]]%in%expr.df$Gene.Name)
    
    mx <- apply(X=expr.df[expr.df$Gene.Name%in%CpGenes[[cp]],tiss.TF], MARGIN=2, 
                FUN=function(tiss.val){
                  
                  # Cp genes not in expr.df + Cp genes in expr.df but have no data for that tissue
                  nND <- nND + sum(is.na(tiss.val))
                  tiss.val <- tiss.val[!is.na(tiss.val)]
                  tiss.val[tiss.val<expr.cutoff] <- 0
                  
                  # Number of Cp genes with data
                  n <- length(tiss.val)
                  if( (nND+n)!=n0 ){
                    print(paste0(cp, ": Not adding up."))
                  }
                  
                  MEAN=mean(tiss.val)
                  SD=sd(tiss.val)
                  SDdivMEAN=SD/MEAN
                  if( !is.finite(SDdivMEAN) ){
                    print(paste0(cp, ": Not finite SDdivMEAN."))
                  }
                  
                  return(
                    
                    c(MEAN=MEAN,
                      MEDIAN=median(tiss.val),
                      SD=SD,
                      # Should not be NaN cause most likely SD would not be 0 in this case
                      SDdivMEAN=SDdivMEAN,
                      # Fraction of Cp genes without data
                      ND=nND/n0,
                      NE=sum(tiss.val==0)/n, 
                      LE=sum(tiss.val>=0.5 & tiss.val<=10)/n,
                      ME=sum(tiss.val>10 & tiss.val<=1000)/n,
                      HE=sum(tiss.val>1000)/n)
                    
                  )
                  
                })
    
    print(paste0(cp, " done!"), quote=FALSE)
    return( cbind.data.frame(Cp=cp, reshape2::melt(mx)) )
    
  })
  
  CPGNEXPR.DF <- do.call("rbind.data.frame", CPGNEXPR.DF)
  colnames(CPGNEXPR.DF) <- c("Cp", "stat", "tissue", "value")
  rownames(CPGNEXPR.DF) <- NULL
  save(CPGNEXPR.DF, file=paste0(out.dir, "/", out.id, ".RData"))
  
} else {
  load(file=paste0(out.dir, "/", out.id, ".RData"))
}

Cp.v <- unique(as.character(CPGNEXPR.DF$Cp))
Cp.v <- Cp.v[Cp.v!="HiC_all"]
Cp.v <- c( "HiC_all", as.character(sort(as.numeric(Cp.v), decreasing=FALSE)) )

CPGNEXPR.DF$Cp <- factor(as.character(CPGNEXPR.DF$Cp), levels=Cp.v)
CPGNEXPR.DF$stat <- as.character(CPGNEXPR.DF$stat)
CPGNEXPR.DF$tissue <- as.character(CPGNEXPR.DF$tissue)

stat.v <- unique(CPGNEXPR.DF$stat)

#-------------------Boxplot

# Boxplot with outliers
pdf(file=paste0(out.dir, "/", out.id, ".pdf"), height=20, width=50)
par(mfrow=c(2,5))

stat1.v <- c("NE", "LE", "ME", "HE", "ND",
             "MEAN", "MEDIAN", "SD", "SDdivMEAN", "ND")

for(stat in stat1.v){
  
  if( stat%in%c("NE", "LE", "ME", "HE", "ND") ){
    ylim <- c(0,1)
  } else {
    ylim <- NULL
  }
  boxplot(value~Cp, outline=TRUE, data=CPGNEXPR.DF[CPGNEXPR.DF$stat==stat,],
          boxwex=0.6, xlab=expression("c"["p"]), ylab=stat, cex.axis=1.5, main="", 
          xaxt="n", col=col.v[stat], ylim=ylim)
  mtext(text=paste0(out.id, "_", length(unique(CPGNEXPR.DF$tissue)), 
                    "tissuesPerBp_", stat), side=3, line=2, cex=1)
  axis(side=1, at=1:length(levels(CPGNEXPR.DF$Cp)), labels=levels(CPGNEXPR.DF$Cp), 
       cex=1.5)
  
}

dev.off()

# Boxplot no outliers
pdf(file=paste0(out.dir, "/", out.id, "_noOutlier.pdf"), height=20, width=50)
par(mfrow=c(2,5))

for(stat in stat1.v){
  
  boxplot(value~Cp, outline=FALSE, data=CPGNEXPR.DF[CPGNEXPR.DF$stat==stat,],
          boxwex=0.6, xlab=expression("c"["p"]), ylab=stat, cex.axis=1.5, main="",
          xaxt="n", col=col.v[stat])
  mtext(text=paste0(out.id, "_", length(unique(CPGNEXPR.DF$tissue)), 
                    "tissuesPerBp_", stat), side=3, line=2, cex=1)
  axis(side=1, at=1:length(levels(CPGNEXPR.DF$Cp)), labels=levels(CPGNEXPR.DF$Cp), 
       cex=1.5)

}

dev.off()

#-------------------Fold change

for( fnx in c("mean", "median") ){
  
  if( is.factor(CPGNEXPR.DF$Cp) ){
    
    eval(parse(text=paste0(
      
      'df <- aggregate(x=CPGNEXPR.DF$value, by=list(CPGNEXPR.DF$stat, CPGNEXPR.DF$Cp), 
                      FUN=', fnx, ')'
      
    )))
    
    
  } else {
    stop("Cp column not factor.")
  }
  
  colnames(df) <- c("stat", "Cp", "value")
 
  # FC of mean of stats relative to Hi-C all
  df$FC <- log2(df$value/df$value[as.character(df$Cp)=="HiC_all"])
  if( sum(is.na(df$FC))!=0 ){
    stop("NA in df$FC")
  }
  
  df$stat <- factor(x=df$stat, levels=stat.v)
  df$Cp <- factor(x=df$Cp, levels=Cp.v)
  p <- ggplot(data=df, aes(x=Cp, y=FC, group=stat, colour=stat)) +
    geom_line(size=2) +
    geom_point(size=2.5) +
    scale_y_continuous(breaks=seq(from=-2, to=0.5, by=0.5), limits=c(-2,0.5)) + 
    scale_colour_manual(values=col.v[levels(df$stat)]) + 
    labs(y=paste0("FC of ", fnx, " of values from tissues"), colour=NULL, 
         title=paste0(out.id, "_reltoHiCall")) +
    bgr2 +
    theme(axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10))
  
  ggsave(filename=paste0(out.dir, "/", out.id, "_aggfnx", fnx, "_line_FC_stat.pdf"), 
         height=10, width=10, plot=p)
  
}

# rm(list=ls()); gc()

  
