################################################################################
# Associate expression (across tissues )of genes co-localised with unique bins 
# per Cp in two ways (see description below.)
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
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
#lib = paste0(home.dir, "/DPhil/lib")
lib = "../lib"
#data.dir = paste0(home.dir, "/Database")
data.dir = "../z_ignore_git/Database"
#wk.dir = paste0(home.dir, "/SahakyanLab/CoreGenomeExplorer")
wk.dir = "z_ignore_git"

exprData.dir = paste0(wk.dir, "/out_cleanExprData")
#cpgenes.dir = paste0(wk.dir, "/txt_CpGenes")
cpgenes.dir = "txt_CpGenes"
out.dir = paste0(wk.dir, "/out_geneExprVsCp")

gene.id = "LTr_ALL" #"ALL"
annofilePath = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
Cp.v = 1:21
gcb = "min2Mb"
expr.cutoff = 0 
src.id = "data2" 
manner = "acrossTissues" # "acrossTissues" # "perTissue"
#statToPlot = c("nND", "nNE", "nLE", "nME", "nHE",
#               "MEAN", "MEDIAN", "SD", "SDdivMEAN")
statToPlot = c("MEAN", "MEDIAN", "SD", "SDdivMEAN")

percNDopenUpperLim = 0.3
plottype = "foldchangeplot" #"boxplot" # "foldchangeplot"
box.outliers = F
dim.boxplot = c(2, 5)
funxforfoldchange = "mean"
plotOnly = T
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/compareTwoDist.R"))
source(paste0(wk.dir, "/lib/plotgeneExprVsCp.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0("expr_", src.id, "_cutoff", expr.cutoff, "_", gene.id, "_", 
                   gcb, "_",  manner)

if(plotOnly==F){
  
  Cp.v <- as.character(sort(Cp.v))
  
  # Get list of unique genes overlapping with all long-range contact regions and 
  # regions per Cp.
  
  cpgenes <- readLines(con=paste0(cpgenes.dir, "/", gcb, "_", gene.id, "_name2"))
  cpheader.ind <- grep(x=cpgenes, pattern=">all_genes_cp_", fixed=T)
  cpgenes <- setNames(object=cpgenes[cpheader.ind + 1L],
                      nm=gsub(x=cpgenes[cpheader.ind], pattern=">all_genes_cp_|_end", 
                              replacement=""))
  
  if( !identical(c("HiC_all", Cp.v), names(cpgenes)) ){
    stop("Cp values are not expected.")
  }
  
  cpgenes <- strsplit(x=cpgenes, split=";", fixed=T)
  
  for( cp in names(cpgenes) ){
    
    if( any(duplicated(cpgenes[[cp]])) ){
      stop(paste0(cp, ": Duplicated genes."))
    }
    
  }
  
  # Expression data
  
  exprDataPath <- paste0(exprData.dir, "/expr_", src.id ,"_cutoff0_", gene.id, ".csv")
  expr.df <- read.csv(file=exprDataPath, header=T, stringsAsFactors=F)
  if( any(duplicated(expr.df$Gene.Name)) ){
    stop("Duplicated gene names in expr.df")
  }
  anno.df <- read.delim(file=annofilePath, header=T, stringsAsFactors=F, sep="\t")[,c("chrom", "name2")]
  expr.df <- expr.df[expr.df$Gene.Name%in%anno.df$name2,]
  rm(anno.df, exprDataPath)
  tiss.TF <- !colnames(expr.df)%in%c("chr", "Gene.Name")
  
  expr.mx <- data.matrix(expr.df[,tiss.TF])
  expr.genes <- expr.df$Gene.Name
  dimnames(expr.mx) <- list(expr.genes, colnames(expr.df)[tiss.TF])
  rm(expr.df)
  
  gc()
  
  # Calculate values per Cp
  
  nTiss <- sum(tiss.TF)
  
  if(manner=="perTissue"){
    margin <- 2
  } else if(manner=="acrossTissues"){
    margin <- 1
  } else {
    stop("Invalid manner argument.")
  }
  
  CPGNEXPR.DF <- sapply(X=names(cpgenes), simplify=F, FUN=function(cp){
    
    nCp <- length(cpgenes[[cp]])
    nCpNDAnytiss <- sum(!cpgenes[[cp]]%in%expr.genes)
    
    mx <- apply(X=expr.mx[expr.genes%in%cpgenes[[cp]],], MARGIN=margin, 
                FUN=function(val.v){
                  
                  nExprND <- sum(is.na(val.v))
                  
                  # Can't be all NA because all genes in expr.df have data in at least 1 tissue.
                  val.v <- val.v[!is.na(val.v)]
                  
                  val.v[val.v < expr.cutoff] <- 0 
                
                  MEAN=mean(val.v)  # Can't be NA 
                  SD=sd(val.v)      # NA with only 1 value
                  SDdivMEAN=SD/MEAN # NA/1 = NA; 1/0 = Inf; 0/0  = NaN 
                  
                  return(
                    
                    c(nExprND=nExprND,
                      MEAN=MEAN,
                      MEDIAN=median(val.v),
                      SD=SD,
                      SDdivMEAN=SDdivMEAN,
                      nNE=sum(val.v < 0.5), 
                      nLE=sum(val.v >= 0.5 & val.v <= 10),
                      nME=sum(val.v > 10 & val.v <= 1000),
                      nHE=sum(val.v > 1000))
                    
                  )
                  
                })
    
    print(paste0(cp, ": mx generated!"), quote=F)
    
    mx <- cbind(nCp=nCp, nCpNDAnytiss=nCpNDAnytiss, t(mx))
    mx <- as.data.frame(mx, stringsAsFactors=F)
    mx <- data.frame(Cp=cp, id=rownames(mx), mx, row.names=NULL)
    
    return(mx)
    
  })
  
  CPGNEXPR.DF <- do.call("rbind.data.frame", CPGNEXPR.DF)
  CPGNEXPR.DF$nTiss <- nTiss
  rownames(CPGNEXPR.DF) <- NULL
  
  if(manner=="perTissue"){
  
    chck.sum <- CPGNEXPR.DF$nCp
    CPGNEXPR.DF$nND <- CPGNEXPR.DF$nCpNDAnytiss + CPGNEXPR.DF$nExprND
  
  } else if(manner=="acrossTissues"){
    
    chck.sum <- CPGNEXPR.DF$nTiss
    CPGNEXPR.DF$nND <- CPGNEXPR.DF$nExprND
    
  }
  
  if( any(rowSums(CPGNEXPR.DF[,c("nND", "nNE", "nLE", "nME", "nHE")])!=chck.sum) ){
    stop("Sums don't add up")
  } else {
    save(CPGNEXPR.DF, file=paste0(out.dir, "/", out.name, ".RData"))
  }
  
  # Checks of possible values
  
  notfin.col <- sapply(X=setdiff(colnames(CPGNEXPR.DF), c("Cp", "id")), 
                       simplify=T, FUN=function(col){
    any(!is.finite(CPGNEXPR.DF[[col]]))
  })
  warning(
    paste0("Non-finite values in: ", 
           paste(names(notfin.col)[notfin.col], collapse=", "))
  )
  
  if( any( !is.finite(CPGNEXPR.DF$SD[!is.na(CPGNEXPR.DF$SD)])) ){
    stop("Unexpected value in SD column.")
  }
  
  
} else {
  
  load(file=paste0(out.dir, "/", out.name, ".RData"))
  
}

# Standardise values
col.std <- ifelse(manner=="perTissue", "nCp", "nTiss")
colToStd <- setdiff(colnames(CPGNEXPR.DF),
                    c("Cp", "id", "nCp", "MEAN", "MEDIAN", "SD", "SDdivMEAN", "nTiss"))
CPGNEXPR.DF[,colToStd] <- CPGNEXPR.DF[,colToStd]/CPGNEXPR.DF[[col.std]]

# 
CPGNEXPR.DF <- CPGNEXPR.DF[CPGNEXPR.DF$nND < percNDopenUpperLim,]

# Plot
plotgeneExprVsCp(CPGNEXPR.DF=CPGNEXPR.DF, 
                 manner=manner, out.dir=out.dir, 
                 out.name=paste0(out.name, "_percNDopenUpperLim", percNDopenUpperLim),
                 plottype=plottype, 
                 plottitle=paste0(out.name, "_percNDopenUpperLim", percNDopenUpperLim),
                 statToPlot=statToPlot, funxforfoldchange=funxforfoldchange,
                 box.outliers=box.outliers, dim.boxplot=dim.boxplot)

# rm(list=ls()); gc()
  
################################################################################
# Associate expression (across tissues )of genes co-localised with unique bins 
# per Cp using two ways.
# The metrics calculated were counts of tissues or genes with no data (nND),
# tissues or genes with no data in the context of expression data only (nExprND),
# not-expressed values (nNE), low-expressed values (nLE), medium-expressed values
# (nME), high-expressed values (nHE), mean, median, sd, and coefficient of 
# variation (sd/mean). Count of no data differs between two ways as described
# below. Metrics unique to each way are described below also. 
# a. Cross-tissue calculation of metrics per gene, such that each boxplot per Cp
# contains gene values calculated across tissues. This looks at the influence of 
# Cp on the expression pattern of each gene across tissues. Calculate metrics
# with values of gene across tissues. In this case, nND and nExprND are the same
# and are equal to the number of tissues without data for that gene. This way
# also needs nTiss, the number of tissues in the expression data, because it
# was used to standardise metrics to convert to fractional values. 
# b. Per-tissue calculation of metrics per Cp gene set, such that each boxplot 
# per Cp contains tissue values calculated for the Cp gene set. The first method 
# looks at the influence of Cp on the expression pattern of each gene across 
# tissues. This investigates the influence of Cp on the collective expression 
# pattern of gene sets per $c_p$, and how this pattern behaves across tissues.
# This method calculates metric per tissue using each cp gene set. Therefore
# it takes note of nCp, the number of unique genes per Cp, and nCpNDAnyTiss,
# the number of Cp genes not in expression data. nExprND in this case is
# equal to number of genes in expression data that have no expression value
# per tissue. nND is then given by nCpNDAnyTiss + nExprND. Standardisation
# to convert to fractional value is done using nCp. 
# Note 1. The code allows to filter values based on the allowable percentage
# of no data, for Method a this filters genes based on the number of tissues
# the gene have a value while for Method b this filters tissue data based
# on the number of cp genes it has data. 
# Note 2. The output for the two methods have the same columns but use
# this description to identify which columns matter to that specific method. 
# Which to present. The original code actually calculated method b likely
# because I was influenced by the hub-based calculation. I later on decided
# to present both because they do measure different things and they differ
# in the coefficient of variation trend, which could be commented on, but
# I realised now that it'd be better to use only method a or the cross-tissue
# calculation because method b assumes that all the genes in a cp gene set
# somehow are coregulated, which likely many not be the case, because
# we're talking about thousands of genes, and the other reason is that
# the fraction of genes without data per cp gene set is around 30-40%,
# which is a lot. Meanwhile the number of tissues without data is <10%
# for all the genes. 

# rm(list=ls()); gc()
