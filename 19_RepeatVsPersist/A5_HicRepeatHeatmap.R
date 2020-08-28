################################################################################
# HicRepeatVsPersist heatmap
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    # directory for intermediate and final files of HiCRepeatExploration
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    # directory for intermediate and final files of HiCRepeatExploration
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "subfam12" # "fam" | "subfam" | "subfam6"
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
PreElmTissDyn.dir = paste0(wk.dir, "/out_HicRepeatHeatmapData/", rep.group)
#hmclustPth = paste0(wk.dir, "/out_HicRepeatHeatmap/hm_famVssubfam_clust.csv")
out.dir = paste0(wk.dir, "/out_HicRepeatHeatmap/", rep.group)
### OTHER SETTINGS #############################################################
# Age rank identifier
out.name = "GiorPubl" 
gcb = "min2Mb"
chr = "chrALL" 
binsize.v = 1L 
# Regenerate ELMTISSDYN.MX?
regenerateData = TRUE
# Lineplot per repeat of average minimum repeat count per Cp
lineplot = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(RColorBrewer)
library(gplots)
library(data.table)
library(ComplexHeatmap)
source(paste0(lib, "/subsetANDbinTableRowCol.R"))

myheatmap <- function(mx=ELMTISSDYN.MX, 
                      colScheme=colorRampPalette(c("darkblue", "blue","lightblue",
                                                   "red", "darkred"))(n = 299),
                      rep.group=rep.group
                      ){
  
  if( grepl(x=rep.group, pattern="subfam", fixed=TRUE) ){
    #heatmap.2(x=mx, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
    #          distfun=function(x) dist(x,method = "euclidean"), scale="none", 
    #          na.rm=FALSE, col=colScheme, na.color="black", trace="none",
    #          margins=c(15,15), cexCol=5, cexRow=0.8, srtCol=0, adjCol=c(0.5,1.5),
    #          key=TRUE, keysize=0.5, key.title=NA, xlab=expression("c"["p"]), 
    #          ylab="Repeat subfamilies in chronological order (oldest to youngest)")
    h1 <- ComplexHeatmap::Heatmap(matrix=mx, col=colScheme, na_col="gray50", 
                                  row_names_gp=gpar(fontsize=2),
                                  cluster_columns=FALSE, 
                                  #row_dend_width=unit(50,"mm")
                                  #, 
                                  cluster_rows=FALSE
                                  )
    return(h1)
  } else if(rep.group=="fam"){
    #heatmap.2(x=mx, Rowv=TRUE, Colv=FALSE, dendrogram="row", 
    #          distfun=function(x) dist(x,method = "euclidean"), scale="none", 
    #          na.rm=FALSE, col=colScheme, na.color="black", trace="none",
    #          margins=c(15,15), cexCol=5, cexRow=2, srtCol=0, adjCol=c(0.5,1.5),
    #          key=TRUE, keysize=0.5, key.title=NA, xlab=expression("c"["p"]), 
    #          ylab="Repeat families")
    h1 <- ComplexHeatmap::Heatmap(matrix=mx, col=colScheme, na_col="gray50",
                                  cluster_columns=FALSE, 
                                  row_dend_width=unit(50,"mm"), 
                                  row_names_gp=gpar(fontsize=15)) 
    return(h1)
  } else {
    stop("Invalid rep.group argument.")
  }
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
coul <- colorRampPalette(c("darkblue", "blue","lightblue","red", "darkred"))(n = 299)
id <- paste0(chr, "_", gcb)
if(regenerateData){
  #---------------------------------------
  # Generate ELMTISSDYN.MX from PREELMTISSDYN.MX
  col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
  agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                      header=TRUE, stringsAsFactors=FALSE)[,col.nme]
  # Load PREELMTISSDYN.MX
  load(file=paste0(PreElmTissDyn.dir, "/", chr, "_", gcb, "_", out.name, 
                   "_PreElmTissDyn.RData"))
  num.contacts.each.ntis <- PREELMTISSDYN.MX[1,]
  ELMTISSDYN.MX <- t( apply(X=PREELMTISSDYN.MX[-1,], MARGIN=1, FUN=function(row){
    #100*(row/num.contacts.each.ntis)
    row/num.contacts.each.ntis
  }) )
  #---------------------------------------
  # Lineplot per repeat of average minimum repeat count per Cp
  ntis <- 1:21
  if(lineplot){
    element <- rownames(ELMTISSDYN.MX)
    for(elm in element){
      vals <- ELMTISSDYN.MX[elm,]
      pdf(file=paste0(out.dir, "/lineplot/", id, "_lineplot_", elm, ".pdf"), 
          width=10, height=10)
      plot(x=ntis, y=vals, col="#55bde6", cex.lab=1.3, cex.axis=1.3, 
           pch=19, lwd=5, xlab="", ylab="", main="")
      lines(x=ntis, y=vals, col="#55bde6")
      # X axis
      mtext(side=1, text=expression("c"["p"]), line=3, cex=1.5)
      # Y axis
      mtext(side=2, text="Contact fraction with non-0 min repeat count", 
            line=2.7, cex=1)
      # Diagram title
      mtext(side=3, text=paste0(id, "_", elm), line=1.5, cex=1.5)
      dev.off()
    } # element for loop end
  }
  pdf(file=paste0(out.dir, "/lineplot/", id, "_densplot_all.pdf"), 
      width=10, height=10)
  d <- density(log10(ELMTISSDYN.MX))
  plot(d, xlab="", ylab="", main="", col="#55bde6")
  polygon(d, col="#55bde6", border="#55bde6")
  # X axis
  mtext(side=1, text=bquote(bold("log"["10"]~"(Contact fraction with non-0 min repeat count)")),
        line=3, cex=1.5)
  # Y axis
  mtext(side=2, text=expression(bold("Density")), line=2.7, cex=1.5)
  # Diagram title
  mtext(side=3, text=paste0(id, "_all_N=", length(rownames(ELMTISSDYN.MX))),
        line=1.5, cex=1.5)
  dev.off()
  #---------------------------------------
  # Normalize values in ELMTISSDYN.MX for making heatmap
  ELMTISSDYN.MX.norm <- ELMTISSDYN.MX
  for(i in 1:dim(ELMTISSDYN.MX)[1]){
    #print(i)
    ELMTISSDYN.MX.norm[i,] <- (ELMTISSDYN.MX[i,] - mean(ELMTISSDYN.MX[i,]))/
      sd(ELMTISSDYN.MX[i,])
  }
  #---------------------------------------
  # Bin normalised values by taking mean; binning only valid for normalised values
  ELMTISSDYN <- sapply(X=binsize.v, simplify=FALSE, FUN=function(binsize){
    if(binsize==1){
      return(ELMTISSDYN.MX.norm)
    } else {
      binned.mx <- subsetANDbinTableRowCol(table=ELMTISSDYN.MX.norm,
                                           toSubsetANDbin=agerank,
                                           defineOrder=TRUE, row=TRUE, 
                                           binBy=binsize)
      return(binned.mx)
    }
  })
  names(ELMTISSDYN) <- paste(paste0(out.name, "binned_NORM"), binsize.v, sep="")
  ELMTISSDYN[[paste0(out.name, "_RAW")]] <- ELMTISSDYN.MX
  rm(ELMTISSDYN.MX, ELMTISSDYN.MX.norm); gc()
  save(ELMTISSDYN, file=paste0(out.dir, "/", chr, "_", gcb, "_ElmTissDyn_",
                               out.name, ".RData"))
} else {
  load(file=paste0(out.dir, "/", chr, "_", gcb, "_ElmTissDyn_",
                   out.name, ".RData"))
}
# Generate heatmap
mx.nme.v <- names(ELMTISSDYN)
for(mx.nme in mx.nme.v){
  mx <- ELMTISSDYN[[mx.nme]]
  pdf(file=paste0(out.dir, "/", id, "_", mx.nme, "_heatmap.pdf"), 
      width=10, height=10)
  h1 <- myheatmap(mx=mx, colScheme=coul, rep.group=rep.group)
  
  #if(rep.group=="fam"){
  #  h2 <- read.csv(file=hmclustPth, header=TRUE)
  #  rownames(h2) <- h2$repName
  #  h2 <- as.matrix(h2[order(h2$rank, decreasing=FALSE),
  #                     c("cluster", "cluster.sf")])
  #  h2[h2=="ClusterI"] <- "inc cluster"
  #  h2[h2=="ClusterII"] <- "dec cluster"
  #  dimnames(h2)[[2]] <- c("fam", "sub-fam")
  #  h2 <- Heatmap(matrix=h2, col=c("gray80", "gray20"), 
  #                row_names_gp=gpar(fontsize=2))
  #  print(h1+h2)
  #} else {
  #   print(h1)
  #}

  print(h1)

  dev.off()
}
print(paste0(chr, ":Heatmap done!"))

# rm(list=ls()); gc()