## FUNCTION ####################################################################
HicRepeatClustering <- function(
  # Chromosome for which the data are to be pulled:
  chr = "chr1",
  # Feature database path and filename suffix:
  featureDBPath = #"/home/alex/Desktop/CHROMSEQ/OUT",
        "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc",
  # Feature database filename suffix
  suffix = "min2Mb",
  # Wheteher to recalculate the ELMTISSDYN.MX again, or use it from the database
  regenerateData = TRUE
){
################################################################################
library(RColorBrewer)
library(gplots)
################################################################################

id <- paste0(chr,"_",suffix)

#########################
if(regenerateData==TRUE){

  load(paste0(featureDBPath,"/",chr,"_MinElm_",suffix,".RData"))

  element.names <- dimnames(MINELM.MX)[[2]][-1]
  hits.len <- length(MINELM.MX[,1])
  unique.ntis <- unique(MINELM.MX[,"ntis"])
  unique.ntis <- unique.ntis[order(unique.ntis)]
  unique.ntis.len <- length(unique.ntis)

  ELMTISSDYN.MX <- matrix(NA, ncol=unique.ntis.len, nrow=length(element.names))
  dimnames(ELMTISSDYN.MX)[[1]] <- element.names
  dimnames(ELMTISSDYN.MX)[[2]] <- unique.ntis

  ###---------------------------
  for(element in element.names){
    print(element, quote=FALSE)
    test1 <- (MINELM.MX[,element]>=1)
    num.contacts.each.ntis <- num.non0cont.each.ntis <- rep(NA,unique.ntis.len)
    k <- 1
    for(ntis in unique.ntis){
      test2 <- (MINELM.MX[,"ntis"]==ntis)
      num.contacts.each.ntis[k] <- sum( test2 )
      num.non0cont.each.ntis[k] <- sum( test2 & test1 )
      k <- k + 1
      rm(test2)
    }; rm(k, test1)

    ELMTISSDYN.MX[element,] <- 100*num.non0cont.each.ntis/num.contacts.each.ntis
  }
  ###---------------------------

  save(ELMTISSDYN.MX, file=paste0(featureDBPath,"/",chr,"_ElmTissDyn_",suffix,".RData"))
}
#########################

if(regenerateData==FALSE){
  load(paste0(featureDBPath,"/",chr,"_ElmTissDyn_",suffix,".RData"))
  unique.ntis <- as.numeric(dimnames(ELMTISSDYN.MX)[[2]])
}

# Row-wise centering and SD-scaling from the ELMTISSDYN.MX matrix
ELMTISSDYN.MX.norm <- ELMTISSDYN.MX
for(i in 1:dim(ELMTISSDYN.MX)[1]){
  #print(i)
  ELMTISSDYN.MX.norm[i,] <- (ELMTISSDYN.MX[i,] - mean(ELMTISSDYN.MX[i,]))/
                                          sd(ELMTISSDYN.MX[i,])
}

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

# ### To plot the lineplot from a single element.
# element <- "MIR"
# quartz()
# #pdf(file=paste0(id,"_lineplot_",element,".pdf"), width=9, height=7)
#   plot(x=unique.ntis, y=ELMTISSDYN.MX[element,],
#        col=rf(max(unique.ntis)), cex=1, pch=21, lwd=3, # ylim=c(0,30)
#        xlab="Num(non-0 tissues/cell lines)",
#        ylab="contact % with >=1 repeat pair",
#        main=paste0(id,"_",element))
#   lines(x=unique.ntis, y=ELMTISSDYN.MX["CR1",])
# #dev.off()

pdf(file=paste0(id,"_CLUSTheatmap.pdf"), width=9, height=9)
# HM <- heatmap.2(ELMTISSDYN.MX, Rowv=TRUE, Colv=FALSE, dendrogram="row",
#                 scale="row", trace="none",  na.color="gray",
#                 xlab="tissue contact freq", ylab="repeats", margins=c(4,10),
#                 #col=colorRampPalette(rf(max(unique.ntis)))(n = 30),
# col=colorRampPalette(c("darkblue", "blue","lightblue","red", "darkred"))(n = 299),
#                 key=TRUE, distfun = function(x) dist(x,method = "euclidean") )
heatmap.2(ELMTISSDYN.MX.norm, Rowv=TRUE, Colv=FALSE, dendrogram="row",
          scale="none", trace="none",  na.color="gray",
          xlab="tissue contact freq", ylab="repeats", margins=c(4,10),
          #col=colorRampPalette(rf(max(unique.ntis)))(n = 30),
col=colorRampPalette(c("darkblue", "blue","lightblue","red", "darkred"))(n = 299),
          key=TRUE, distfun = function(x) dist(x,method = "euclidean") )
dev.off()

pdf(file=paste0(id,"_CLUSTheatmapRAW.pdf"), width=9, height=9)
heatmap.2(ELMTISSDYN.MX, Rowv=TRUE, Colv=FALSE, dendrogram="row",
          scale="none", trace="none",  na.color="gray",
          xlab="tissue contact freq", ylab="repeats", margins=c(4,10),
          #col=colorRampPalette(rf(max(unique.ntis)))(n = 30),
col=colorRampPalette(c("darkblue", "blue","lightblue","red", "darkred"))(n = 299),
# col=rf(100), #col=colorRampPalette(c("magenta", "green"))(n = 299), # col=redgreen,
          key=TRUE, distfun = function(x) dist(x,method = "euclidean") )
dev.off()


set.seed(845)
# centers=7 is the best from silhouette analyses.
# The above analysis was done via the command:
# fviz_nbclust(ELMTISSDYN.MX.norm, kmeans, method = "silhouette")
clust <- kmeans(x=ELMTISSDYN.MX.norm, centers=7, iter.max = 1000, nstart = 100, trace=FALSE)
pdf(file=paste0(id,"_CLUSTlineplot.pdf"), width=12, height=12)
par(mfrow=c(3,3))
for(cl in 1:7){
  elements.inclust <- names(which(clust$cluster==cl))
  y.range <- range(ELMTISSDYN.MX.norm[elements.inclust,])

  plot(NA, ylim=y.range, xlim=range(unique.ntis),
       col=rf(max(unique.ntis)), cex=1, pch=21, lwd=3, # ylim=c(0,30)
       xlab="Num(non-0 tissues/cell lines)",
       ylab="Scaled contact % with >=1 repeat pair",
       cex.main=0.4,
       main=paste0(id,"_clust",cl,"_",length(elements.inclust),"\n",paste0(elements.inclust,collapse=";")) )
  for(i in 1:21){abline(v=i, col="grey", lty="dotted")}
  for(elm in elements.inclust){
    lines(x=unique.ntis, y=ELMTISSDYN.MX.norm[elm,], col="grey", lwd=3)
  }
  lines(x=unique.ntis, y=clust$centers[cl,], col="navy", lwd=4)
}
dev.off()

print("HicRepeatClustering is DONE!", quote=FALSE)

}
################################################################################
suppressWarnings(suppressPackageStartupMessages(library(compiler)))
HicRepeatClustering <- cmpfun(HicRepeatClustering, options=list(suppressUndefined=TRUE))


# Execution of the function written above...
################################################################################
HicRepeatClustering(
  # Chromosome for which the data are to be pulled:
  chr = "chr1",
  # Feature database path and filename suffix:
  featureDBPath = #"/home/alex/Desktop/CHROMSEQ/OUT",
        "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc",
  # Feature database filename suffix
  suffix = "min2Mb",
  # Wheteher to recalculate the ELMTISSDYN.MX again, or use it from the database
  regenerateData = FALSE
)




#library(ggplot2)
#p <- ggplot(df, aes(ntis,valsum)) +
#      stat_bin2d(binwidth=c(1,10)) +
#      scale_fill_gradientn(colours=r)

#> names(PERSIST.MX)
#[1] "hits"    "ntis"    "valsum"  "control"

#> PERSIST.MX$hits[1,]
#       i  j Co Hi Lu LV RV Ao PM Pa Sp Li SB AG Ov Bl MesC MSC NPC TLC ESC FC LC
#179619 2 54  0  0  0  0  0  0  0  0  0  0  0  0  0  0    0   1   0   0   0  1  0

#> PERSIST.MX$ntis
#[1]  2  7  2 10 14  7 12 ...

#> PERSIST.MX$valsum
#[1]   2  14   2  23  35  ...

#> PERSIST.MX$control[1,]
#i  j value
#172840 1 52     0


#> BINKMER4.MX[1,]
#bins startpos   endpos     AAAA     AAAC     AAAG     AAAT
#2    40001    80000       NA       NA       NA       NA   # NA means there is at least one N

#> BINREP.MX[1,]
#bins  startpos  endpos    ERV1    L1  ...
#2          40001          80000   0   ...

#> names(BINREP.MX[1,])
#[1] "bins"           "startpos"       "endpos"         "ERV1"           "L1"
#[6] "Alu"            "ERVL"           "hAT-Charlie"    "ERVL-MaLR"      "MIR"
#[11] "Satellite"      "ERVK"           "Low_complexity" "L2"             "TcMar-Tigger"
#[16] "Simple_repeat"  "hAT-Blackjack"  "Unknown"        "centr"          "RTE"
#[21] "TcMar-Mariner"  "Gypsy?"         "CR1"            "LTR"            "hAT-Tip100"
#[26] "hAT"            "Gypsy"          "rRNA"           "DNA"            "tRNA"
#[31] "TcMar?"         "srpRNA"         "TcMar-Tc2"      "ERV"            "snRNA"
#[36] "Helitron"       "hAT?"           "RTE-BovB"       "SINE"           "ERVL?"
#[41] "PiggyBac"       "LTR?"           "DNA?"           "Other"          "Merlin"
#[46] "MuDR"           "SINE?"          "Deu"            "RNA"            "scRNA"
#[51] "TcMar"          "Dong-R4"        "PiggyBac?"      "Unknown?"       "Penelope?"
#[56] "L1?"            "Helitron?"      "telo"

#
# [9] "chr10_BinKmer4_min05Mb.RData" "chr10_BinKmer4_min2Mb.RData"
# [11] "chr10_BinKmer7_min05Mb.RData" "chr10_BinKmer7_min2Mb.RData"
# [13] "chr10_BinRep_min05Mb.RData"   "chr10_BinRep_min2Mb.RData"
# [15] "chr10_Persist_min05Mb.RData"  "chr10_Persist_min2Mb.RData"
#
