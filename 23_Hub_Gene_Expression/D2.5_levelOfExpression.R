################################################################################
# Categorise hubs based on level of expression of genes.  EMBL-EBI defines the 
# categories as not expressed (TPM/FPKM<0.5), lowly-expressed (TPM/FPKM within 
# [0.5,10]), mediumly-expressed (TPM/FPKM within (10,1000]) and highly-expressed 
# (TPM/FPKM>1000). 
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
    data.dir = "/Users/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
hub.dir = paste0(wk.dir, "/out_hubfile/All_LTr/acceptable_n")
out.dir = paste0(wk.dir, "/out_levelOfExpression/acchubs/All_LTr")

#hub.dir = paste0(wk.dir, "/out_hubfile/acceptable_n")
#out.dir = paste0(wk.dir, "/out_levelOfExpression/acchubs")
### OTHER SETTINGS #############################################################
# Protein-coding gene expression data (bulk RNA-seq)
exprData.dir = paste0(wk.dir, "/out_cleanExprData")
exprData.suffix = "cutoff0_LTr_ALL"
src.id.v = c("data1", "data2")
gcb = "min2Mb"
hub.id = "All_topCP3_gapBin50"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(gcb, "_", hub.id)

EXPR.DF <- list()
for(src.id in src.id.v){
  
  HUB <- list.files(path=hub.dir, pattern=hub.id, recursive=F)
  HUB <- HUB[grepl(x=HUB, pattern=gcb)]
  if(length(HUB)==0){
    stop( paste0(src.id, ": No hubs.") ) 
  }
  
  expr.df <- read.csv(file=paste0(exprData.dir, "/expr_", src.id, "_", exprData.suffix, 
                                  ".csv"), header=T, stringsAsFactors=F)
  tiss.v <- sort(colnames(expr.df[,!colnames(expr.df)%in%c("chr", "Gene.Name")]))
  
  for(hub in HUB){
    
    hub <- gsub(x=hub, pattern=".csv", replacement="", fixed=T)
    
    genes <- read.csv(file=paste0(hub.dir, "/", hub, ".csv"), header=T, 
                      stringsAsFactors=F)[,"gene"]
    if( all(is.na(genes)) ){
      print(paste0(src.id, ": No genes for ", hub, "."), quote=FALSE)
      next
    }
    genes <- unlist(strsplit(x=genes, split=";"))
    genes <- unique(genes[!is.na(genes)])
    
    id <- paste0(src.id, "_", hub)
    EXPR.DF[[id]] <- cbind.data.frame(src.id=src.id, 
                                      hub=gsub(hub, pattern=paste0(gcb, "_|", hub.id, "_"), replacement=""),
                                      reshape2::melt(expr.df[expr.df$Gene.Name%in%genes,c("Gene.Name", tiss.v)],
                                                     id="Gene.Name"),
                                      stringsAsFactors=F)
    
    if( !identical(levels(EXPR.DF[[id]]$variable), tiss.v) ){
      stop(paste0(id, ": Error in melting table."))
    }
    
    print(paste0(id, " done!"), quote=F)
    
  } # HUB for loop end
  
} # src.id.v for loop end

EXPR.DF <- do.call("rbind", EXPR.DF)
rownames(EXPR.DF) <- NULL
EXPR.DF$value <- log10(EXPR.DF$value)
EXPR.DF$hub <- factor(x=as.character(EXPR.DF$hub), 
                      levels=sort(unique(as.character(EXPR.DF$hub)), decreasing=F))
EXPR.DF$src.id <- factor(x=as.character(EXPR.DF$src.id), 
                         levels=sort(src.id.v, decreasing=F))

HUB.len <- length(levels(EXPR.DF$hub))

EXPR.DF$hub <- as.character(EXPR.DF$hub)
EXPR.DF$hub <- gsub(x=EXPR.DF$hub, pattern="_", replacement="\n", fixed=TRUE)
EXPR.DF$hub <- factor(x=as.character(EXPR.DF$hub), 
                      levels=sort(unique(as.character(EXPR.DF$hub)), decreasing=F))

# Boxplot

col.v <- c("gray80", "gray50") #c("aquamarine3", "coral") #c("white", "gray")
leg.id <- paste(paste(levels(EXPR.DF$src.id), col.v, sep="-"), collapse=",")

pdf(file=paste0(out.dir, "/", out.id, "_levExpr_bp.pdf"), 
    width=40, height=10)
par(mar = c(6.1, 6.1, 4.1, 4.1))

boxplot(value~src.id*hub, outline=T, data=EXPR.DF, boxwex=0.6, xlab=NULL, 
        ylab=NULL, cex.axis=2.5, col=col.v, xaxt="n", yaxt="n", 
        main=paste0(out.id, "_legend:", leg.id, "_log10(expression value) all tissues")
        ) 
axis(side=1, at=seq(1.5, HUB.len*2, 2), labels=levels(EXPR.DF$hub), cex.axis=3, 
     mgp=c(8, 4.5, 0))
axis(side=2, las=2, cex.axis=3)
abline(v=seq(0.5, HUB.len*2+0.5, by=2), lty=1, lwd=2, col="gray50") 
abline(h=log10(c(0.5, 10, 1000)), lty=2, lwd=5, col=c("#2171B5", "#9ECAE1", "#A50F15"))
legend(x="topright", legend=levels(EXPR.DF$src.id), col=col.v, lty=1, lwd=5, cex=1, bty="n")

dev.off()

# rm(list=ls()); gc()

