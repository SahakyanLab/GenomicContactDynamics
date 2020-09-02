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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    # directory for intermediate and final files of HiCRepeatExploration
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "fam" # "fam" | "subfam" | "subfam6"
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
PreElmTissDyn.dir = paste0(wk.dir, "/out_HicRepeatHeatmapData/", rep.group)
# hmclustPth = paste0(wk.dir, "/out_HicRepeatHeatmap/hm_famVssubfam_clust.csv")
out.dir = paste0(wk.dir, "/out_HicRepeatHeatmap")
### OTHER SETTINGS #############################################################
# Age rank identifier
out.name = "GiorPubl" 
gcb = "min2Mb"
chr = "chrALL" 
# Regenerate ELMTISSDYN?
regenerateData = TRUE
# Lineplot per repeat of average minimum repeat count per Cp
lineplot = TRUE
cluster.TF = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(viridis)
library(gplots)
library(data.table)
library(ComplexHeatmap)


myheatmap <- function(mx=ELMTISSDYN.MX, 
                      colScheme=coul <- viridis::viridis(n=299),
                      rep.group=rep.group, mx.nme="raw", cluster=FALSE
                      ){
  if(mx.nme=="raw"){
    at.v=seq(from=0, to=1, by=0.2)
  } 
  
  if( grepl(x=rep.group, pattern="subfam", fixed=TRUE) ){
    h1 <- ComplexHeatmap::Heatmap(matrix=mx, col=colScheme, na_col="gray50", 
                                  cluster_columns=FALSE, cluster_rows=cluster,
                                  row_names_gp=gpar(fontsize=2),
                                  heatmap_legend_param=list(title=mx.nme))
    if(mx.nme=="raw"){
      h1 <- ComplexHeatmap::Heatmap(matrix=mx, col=colScheme, na_col="gray50", 
                                    cluster_columns=FALSE, cluster_rows=cluster,
                                    row_names_gp=gpar(fontsize=2),
                                    heatmap_legend_param=list(title=mx.nme, at=at.v))
    }
    return(h1)
  } else if(rep.group=="fam"){
    
    h1 <- ComplexHeatmap::Heatmap(matrix=mx, col=colScheme, na_col="gray50",
                                  cluster_columns=FALSE, cluster_rows=cluster, 
                                  row_dend_width=unit(50,"mm"), 
                                  row_names_gp=gpar(fontsize=15),
                                  heatmap_legend_param=list(title=mx.nme))
    if(mx.nme=="raw"){
      h1 <- ComplexHeatmap::Heatmap(matrix=mx, col=colScheme, na_col="gray50",
                                    cluster_columns=FALSE, cluster_rows=cluster, 
                                    row_dend_width=unit(50,"mm"), 
                                    row_names_gp=gpar(fontsize=15),
                                    heatmap_legend_param=list(title=mx.nme, at=at.v))
    } 
    return(h1)
  } else {
    stop("Invalid rep.group argument.")
  }
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * 
################################################################################
out.dir <- paste0(out.dir, "/", rep.group)
if( !dir.exists(out.dir) ){
  dir.create(out.dir)
}
linep.dir <- paste0(out.dir, "/lineplot")
if( !dir.exists(linep.dir) ){ dir.create(linep.dir) }

coul <- viridis::viridis(n=100)

id <- paste0(chr, "_", gcb)
ELMTISSDYN <- list()
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
  ELMTISSDYN[["raw"]] <- t( apply(X=PREELMTISSDYN.MX[-1,], MARGIN=1, FUN=function(row){
    # 100*(row/num.contacts.each.ntis)
    row/num.contacts.each.ntis
  }) )
  #---------------------------------------
  # Lineplot per repeat of average minimum repeat count per Cp
  ntis <- 1:21
  if(lineplot){
    element <- rownames(ELMTISSDYN$raw)
    for(elm in element){
      vals <- ELMTISSDYN$raw[elm,]
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
  pdf(file=paste0(linep.dir,  "/", id, "_densplot_all.pdf"), 
      width=10, height=10)
  d <- density(log10(ELMTISSDYN$raw))
  plot(d, xlab="", ylab="", main="", col="#55bde6")
  polygon(d, col="#55bde6", border="#55bde6")
  # X axis
  mtext(side=1, text=bquote(bold("log"["10"]~"(Contact fraction with non-0 min repeat count)")),
        line=3, cex=1.5)
  # Y axis
  mtext(side=2, text=expression(bold("Density")), line=2.7, cex=1.5)
  # Diagram title
  mtext(side=3, text=paste0(id, "_all_N=", length(rownames(ELMTISSDYN$raw))),
        line=1.5, cex=1.5)
  dev.off()
  #---------------------------------------
  # Normalize (z-score) values in ELMTISSDYN$raw for making heatmap
  ELMTISSDYN[["norm"]] <- ELMTISSDYN$raw
  for(i in 1:dim(ELMTISSDYN$raw)[1]){
    ELMTISSDYN$norm[i,] <- (ELMTISSDYN$raw[i,] - mean(ELMTISSDYN$raw[i,]))/
      sd(ELMTISSDYN$raw[i,])
  }
  ELMTISSDYN[["fc"]] <- log2(ELMTISSDYN$raw/ELMTISSDYN$raw[,"1"])
  save(ELMTISSDYN, file=paste0(out.dir, "/", chr, "_", gcb, "_ElmTissDyn_",
                               out.name, ".RData"))
} else {
  load(file=paste0(out.dir, "/", chr, "_", gcb, "_ElmTissDyn_", out.name, 
                   ".RData"))
}

# Generate heatmap
mx.nme.v <- names(ELMTISSDYN)
for(mx.nme in mx.nme.v){
  mx <- ELMTISSDYN[[mx.nme]]
  
  pdf(file=paste0(out.dir, "/", id, "_", mx.nme, "_heatmap.pdf"), 
      width=10, height=10)
  
  drp.rw <- apply( X=mx, MARGIN=1, FUN=function(rw) any(!is.finite(rw)) )
  mx[!is.finite(mx)] <- NA
  
  if(cluster.TF){
    # Rows with non-finite values; exclude in clustering then append later to
    # clustered heatmap.
    h <- myheatmap(mx=mx[!drp.rw,], colScheme=coul, rep.group=rep.group, 
                   mx.nme=mx.nme, cluster=cluster.TF)
    h2 <- myheatmap(mx=mx[drp.rw,], colScheme=coul, rep.group=rep.group, 
                    mx.nme=mx.nme, cluster=FALSE)
    h <- h %v% h2
    
  } else {
    h <- myheatmap(mx=mx, colScheme=coul, rep.group=rep.group, 
                   mx.nme=mx.nme, cluster=cluster.TF)
  }
 
  print(h)
  
  dev.off()
}

# rm(list=ls()); gc()

