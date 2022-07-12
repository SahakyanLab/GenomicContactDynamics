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
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/CoreGenomeExplorer")

exprData.dir = paste0(wk.dir, "/out_cleanExprData")
pairCor.dir = paste0(wk.dir, "/out_coexpression_pairCor")
pairCp.dir = paste0(wk.dir, "/out_coexpression_pairCp")
pairHub.dir = paste0(wk.dir, "/out_coexpression_pairHub")
out.dir = paste0(wk.dir, "/z_ignore_git/out_coexpression_vsCp")

gene.id = "LTr_ALL" #"ALL"
anno.file = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
bin.len = 40000
chr.v = paste0("chr", c(1:22, "X"))
cor.meth = "pearson"
expr.cutoff = 0.5
gcb = "min2Mb"
genepairMaxCp.v = NULL # If NULL, allow all Cp values
genepairMaxCp.id = "1To21" 
hub.id = "min2Mb_All_topCP3_gapBin50"
percNDclosedUpperLim = 0.5
plotOnly = T
src.id = "data2"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggsci)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/generateCORCPDF.R"))
source(paste0(wk.dir, "/lib/compareTwoDist.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0("expr_", src.id, "_cutoff", expr.cutoff, "_", gene.id, 
                   "_", gcb, "_binlen", bin.len, "_", cor.meth)

if(plotOnly==F){
  
  generateCORCPDF(chr.v=chr.v, pairCor.dir=pairCor.dir, src.id=src.id, 
                  expr.cutoff=expr.cutoff, gene.id=gene.id, cor.meth=cor.meth, 
                  pairCp.dir=pairCp.dir, gcb=gcb, bin.len=bin.len, 
                  pairHub.dir=pairHub.dir, hub.id=hub.id,
                  out.dir=out.dir, out.name=out.name)
                              
  
} else {
  load(paste0(out.dir, "/", out.name, "_plot.RData"))
}

out.name <- paste0(out.name, "_percNDclosedUpperLim", percNDclosedUpperLim,
                   "_genepairMaxCp", genepairMaxCp.id)

#
CORCP.DF <- CORCP.DF[CORCP.DF$frTisspairAtleast1NA <= percNDclosedUpperLim,]
CORCP.DF$frTisspairAtleast1NA <- NULL

#
CORCP.DF <- CORCP.DF[!is.na(CORCP.DF$corval),] 
if( !is.null(genepairMaxCp.v) ){
  CORCP.DF <- CORCP.DF[CORCP.DF$maxCp%in%genepairMaxCp.v,]
}

# Correlation coefficient vs. Cp

Cp.v <- as.character(sort(unique(CORCP.DF$maxCp)))
CORCP.DF$maxCp <- factor(x=as.character(CORCP.DF$maxCp), levels=Cp.v)

png(filename=paste0(out.dir, "/", out.name, "_VsMaxCp.png"), 
    height=300*10, width=300*10, res=300)

boxplot(formula=corval~maxCp, data=CORCP.DF[!is.na(CORCP.DF$maxCp),], outline=F, xlab="Cp", 
        ylab=cor.meth, main=out.name, col="#FDC776", cex.main=0.8)

dev.off()

write.csv(x=stack(table(as.numeric(as.character(CORCP.DF$maxCp)), useNA="always")), 
          row.names=F, file=paste0(out.dir, "/", out.name, "_NdatapointsPerCp.csv"))

CORCP.DF$maxCp <- as.character(CORCP.DF$maxCp)
TEST <- compareTwoDist(x=CORCP.DF$corval[CORCP.DF$maxCp=="1"],
                       y=CORCP.DF$corval[CORCP.DF$maxCp=="21"])

# Correlation coefficient non-hub vs. hub

CORCP.DF <- CORCP.DF[!is.na(CORCP.DF$hubTF),]
CORCP.DF$hubTF <- as.character(CORCP.DF$hubTF)

TEST <- compareTwoDist(x=CORCP.DF$corval[CORCP.DF$hubTF=="0"],
                       y=CORCP.DF$corval[CORCP.DF$hubTF=="1"])

mean.df <- aggregate(x=CORCP.DF$corval, by=list(CORCP.DF$hubTF), FUN=mean)
rownames(mean.df) <- mean.df$Group.1
median.df <- aggregate(x=CORCP.DF$corval, by=list(CORCP.DF$hubTF), FUN=median)
rownames(median.df) <- median.df$Group.1

test.id <- paste0("\nmean_0;1=", paste(mean.df[c("0","1"),"x"], collapse=";"),
                  "_median_0;1=", paste(median.df[c("0","1"),"x"], collapse=";"),
                  "\ntwosided", 
                  "_ttest_", TEST$t$p.value, 
                  "_mwtest_", TEST$mw$p.value,
                  "_kstest_", TEST$ks$p.value)

col.v = c(`0`="gray60", `1`="tomato")
count01.v <- table(CORCP.DF$hubTF)

##

backup <- CORCP.DF
CORCP.DF <- CORCP.DF[CORCP.DF$maxCp%in%as.character(c(1:3, 19:21)),]
CORCP.DF$maxCp <- as.numeric(CORCP.DF$maxCp)
CORCP.DF$maxCp[CORCP.DF$maxCp <=3] <- 1
CORCP.DF$maxCp[CORCP.DF$maxCp >=19] <- 21
CORCP.DF$maxCp <- as.character(CORCP.DF$maxCp)
CORCP.DF$id <- paste0(CORCP.DF$maxCp, "_", CORCP.DF$hubTF)

##
p <- ggplot(data=CORCP.DF, aes(x=corval)) +  #, group=hubTF)) +
  geom_density(aes(colour=id)) +
  scale_y_continuous(limits=c(0,2)) + 
  #scale_colour_manual(values=col.v) +
  #scale_fill_manual(values=col.v) + 
  labs(title=paste0(out.name, 
                    "_N0=", count01.v["0"],
                    "_N1=", count01.v["1"], test.id), 
       #x=paste0("(", cor.meth, "-mean)/sd")
       x=cor.meth) +
  bgr2

boxplot(formula=corval~id, data=CORCP.DF[!is.na(CORCP.DF$maxCp),], outline=F, xlab="Cp", 
        ylab=cor.meth, main=out.name, col="#FDC776", cex.main=0.8)
  
ggsave(filename=paste0(out.dir, "/", out.name, "_nonhub0VsHub1.png"), 
       height=10, width=10, plot=p)

# Add gene expression distribution for non-hub and hub gene pairs

g.ind.nhub <- unique(unlist( CORCP.DF[CORCP.DF$hubTF=="0",c("g1.ind", "g2.ind")] ))
g.ind.hub <- unique(unlist( CORCP.DF[CORCP.DF$hubTF=="1",c("g1.ind", "g2.ind")] ))

anno.v <- read.delim(file=anno.file, header=T, stringsAsFactors=F, sep="\t")[,"name2"]
exprData.file <- paste0(exprData.dir, "/expr_", src.id ,"_cutoff0_", gene.id, ".csv")
expr.df <- read.csv(file=exprData.file, header=T, stringsAsFactors=F, row.names="Gene.Name")
expr.df <- expr.df[,colnames(expr.df)!="chr"]

g.expr.nhub <- unlist(expr.df[anno.v[g.ind.nhub],], use.names=F)
g.expr.hub <- unlist(expr.df[anno.v[g.ind.hub],], use.names=F)

bp.df <- data.frame(hubTF=c( rep("0", times=length(g.expr.nhub)), 
                             rep("1", times=length(g.expr.hub)) ),
                    expression=c(g.expr.nhub, g.expr.hub),
                    stringsAsFactors=T)

bp.df$expression[bp.df$expression < expr.cutoff] <- 0

TEST <- compareTwoDist(x=bp.df$expression[bp.df$hubTF=="0"],
                       y=bp.df$expression[bp.df$hubTF=="1"])

mean.df <- aggregate(x=bp.df$expression, by=list(bp.df$hubTF), FUN=mean, na.rm=T)
rownames(mean.df) <- mean.df$Group.1
median.df <- aggregate(x=bp.df$expression, by=list(bp.df$hubTF), FUN=median, na.rm=T)
rownames(median.df) <- median.df$Group.1

test.id <- paste0("\nmean_0;1=", paste(mean.df[c("0","1"),"x"], collapse=";"),
                  "_median_0;1=", paste(median.df[c("0","1"),"x"], collapse=";"),
                  "\ntwosided", 
                  "_ttest_", TEST$t$p.value, 
                  "_mwtest_", TEST$mw$p.value,
                  "_kstest_", TEST$ks$p.value)

count01.v <- table(CORCP.DF$hubTF)

png(filename=paste0(out.dir, "/", out.name, "_bpExpression_nonhub0VsHub1.png"), 
    height=300*10, width=300*10, res=300)

boxplot(formula=expression~hubTF, data=bp.df, outline=F, xlab="Non-hub0; Hub1", 
        ylab="expression in TPM; all tissues", ylim=c(0,60),
        col=adjustcolor(col.v[levels(bp.df$hubTF)], alpha.f=0.5), cex.main=0.8,
        main=paste0(out.name, "_N0=", count01.v["0"], "_N1=", count01.v["1"], test.id))
        
dev.off()

# rm(list=ls()); gc()
