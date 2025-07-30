################################################################################
# Plot gene pair correlation vs. Cp and vs. presence in hub. For correlation vs.
# max Cp, entries with missing max Cp and/or correlation coefficient are not
# counted and displayed. For correlation vs. hub, entries with missing max Cp
# are considered because genes may not be connected by a contact, they can still
# be on the same hub.
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
#lib = paste0(home.dir, "/DPhil/lib")
#data.dir = paste0(home.dir, "/Database")
#wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/23_Hub_Gene_Expression")
lib = "../lib"
data.dir = wk.dir = "z_ignore_git"

exprData.dir = paste0(wk.dir, "/out_cleanExprData")
pairCor.dir = paste0(wk.dir, "/out_coexpression_pairCor")
pairCp.dir = paste0(wk.dir, "/out_coexpression_pairCp")
pairHub.dir = paste0(wk.dir, "/out_coexpression_pairHub")
out.dir = paste0(wk.dir, "/out_coexpression_vsCp")

gene.id = "LTr_ALL" #"ALL"
anno.file = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
bin.len = 40000
chr.v = paste0("chr", c(1:22, "X"))
cor.meth = "pearson" #"spearman"
expr.cutoff = 0 #0.5
gcb = "min2Mb"
genepairMaxCp.v = NULL # If NULL, allow all Cp values
genepairMaxCp.id = "1To21" 
hub.id = "min2Mb_All_topCP3_gapBin50"
percNDclosedUpperLim = 0.5 #0.5
plotOnly = F
src.id = "data2"

# For correlation vs. group plots
Cp.v = 1:21
dynamic.Cp = 1:3
persistent.Cp = 19:21
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggsci)
library(RColorBrewer)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/makeViolinPlot.R"))
source(paste0(lib, "/compareManyDist.R"))
source(paste0(lib, "/compareTwoDist.R"))
source(paste0(wk.dir, "/lib/vsCp_generateCORCPDF.R"))
source(paste0(wk.dir, "/lib/vsCp_makecorVsCpPlot.R"))
source(paste0(wk.dir, "/lib/vsCp_makecorVsGroupPlot.R"))
source(paste0(wk.dir, "/lib/vsCp_makeGeneExprVsGroupPlot.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0("expr_", src.id, "_cutoff", expr.cutoff, "_", gene.id, 
                   "_", gcb, "_binlen", bin.len, "_", cor.meth)

if(plotOnly==F){
  
  CORCP.DF <- generateCORCPDF(chr.v=chr.v, pairCor.dir=pairCor.dir, src.id=src.id, 
                  expr.cutoff=expr.cutoff, gene.id=gene.id, cor.meth=cor.meth, 
                  pairCp.dir=pairCp.dir, gcb=gcb, bin.len=bin.len, 
                  pairHub.dir=pairHub.dir, hub.id=hub.id,
                  out.dir=out.dir, out.name=out.name)
                            
} else {
  load(paste0(out.dir, "/", out.name, "_plot.RData"))
}

# Record pertinent numbers
paircount.mx <- matrix(data=NA, nrow=6, ncol=length(Cp.v), 
                       dimnames=list(c("WithMaxCp", "WithDataMin1Tiss", 
                                       paste0("WithLEQ", percNDclosedUpperLim, "percNDclosedUpperLim"),
                                       "WithFiniteCor", "WithSigCorBHless0.05", "SigCorWithHubTF"), Cp.v))

# Possible gene pairs (with or without maxCp)
possible.pair.count <- nrow(CORCP.DF)

CORCP.DF <- CORCP.DF[!is.na(CORCP.DF$maxCp), ]
# Pairs with maxCp
for(Cp in Cp.v){
  paircount.mx[1,as.character(Cp)] <- sum(CORCP.DF$maxCp == Cp)
}

CORCP.DF <- CORCP.DF[CORCP.DF$frTisspairAtleast1NA < 1, ]
# Pairs with data in at least 1 tissue
for(Cp in Cp.v){
  paircount.mx[2,as.character(Cp)] <- sum(CORCP.DF$maxCp == Cp)
}

# Take only correlations with tissue data <= percNDclosedUpperLim 

CORCP.DF <- CORCP.DF[CORCP.DF$frTisspairAtleast1NA <= percNDclosedUpperLim,]
#CORCP.DF$frTisspairAtleast1NA <- NULL
# Pairs satisfying percNDclosedUpperLim
for(Cp in Cp.v){
  paircount.mx[3,as.character(Cp)] <- sum(CORCP.DF$maxCp == Cp)
}

#

CORCP.DF <- CORCP.DF[is.finite(CORCP.DF$corval), ]
sum(!is.finite(CORCP.DF$corpval)) # Pairs with non-finite correlation pvalues

# Pairs with finite correlation
for(Cp in Cp.v){
  paircount.mx[4,as.character(Cp)] <- sum(CORCP.DF$maxCp == Cp)
}

##
#CORCP.DF <- CORCP.DF[!is.na(CORCP.DF$corval),] 

# Adjust p-values

CORCP.DF$corpval <- p.adjust(CORCP.DF$corpval, method="BH")
CORCP.DF <- CORCP.DF[CORCP.DF$corpval < 0.05, ]
# Pairs with significant correlations
for(Cp in Cp.v){
  paircount.mx[5,as.character(Cp)] <- sum(CORCP.DF$maxCp == Cp)
}

#boxplot(frTisspairAtleast1NA~maxCp, data=tmp, outline=F)
#

# count.df <- table(CORCP.DF[,c("maxCp", "hubTF")], useNA="always")
# colnames(count.df) <- paste0("hubTF", colnames(count.df))
# rownames(count.df) <- paste0("maxCp", rownames(count.df))
# write.csv(x=count.df, row.names=T, 
#           file=paste0(out.dir, "/", out.name, "_percNDclosedUpperLim", 
#                       percNDclosedUpperLim, "_NAcorvalueRemoved_NdatapointsPerCp.csv"))

#
out.name <- paste0(out.name, "_percNDclosedUpperLim", percNDclosedUpperLim,
                   "_genepairMaxCp", genepairMaxCp.id)

#CORCP.DF <- CORCP.DF[!is.na(CORCP.DF$maxCp),]
if( !is.null(genepairMaxCp.v) ){
  CORCP.DF <- CORCP.DF[CORCP.DF$maxCp%in%genepairMaxCp.v,]
}

# Correlation coefficient vs. Cp

makecorVsCpPlot(CORCP.DF=CORCP.DF, CpsToCompare=list(1:3, 19:21), out.dir=out.dir, out.name=out.name)

# Correlation coefficient distributions of dynamic/persistent max Cp and non-hub/hub gene pairs

out.id <- paste0(paste(range(dynamic.Cp), collapse="-"), "_",
                 paste(range(persistent.Cp), collapse="-"))

CORCP.DF <- CORCP.DF[!is.na(CORCP.DF$hubTF),]
# Pairs with significant correlations and nonhub, hub information
for(Cp in Cp.v){
  paircount.mx[6,as.character(Cp)] <- sum(CORCP.DF$maxCp == Cp)
}

#
write.csv(paircount.mx, file=paste0(out.dir, "/", out.name, "_paircountmx", "_countPossibleGenePair", 
                                    possible.pair.count, ".csv"), row.names=T)

CORCP.DF <- makecorVsGroupPlot(CORCP.DF=CORCP.DF, dynamic.Cp=dynamic.Cp, 
                               persistent.Cp=persistent.Cp, 
                               out.dir=out.dir, out.name=paste0(out.name, "_", out.id))

# Gene expression distribution per gene pair group

makeGeneExprVsGroupPlot(CORCP.DF=CORCP.DF, anno.file=anno.file, 
                        exprData.dir=exprData.dir, src.id=src.id, 
                        gene.id=gene.id, expr.cutoff=expr.cutoff, 
                        out.dir=out.dir, out.name=paste0(out.name, "_", out.id),
                        removeOutliers = F, addmean = F)

# rm(list=ls()); gc()
