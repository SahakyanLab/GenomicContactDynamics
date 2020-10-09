################################################################################
# Heatmap summarising VUS of subj-ref comparisons
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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
VUScsvPath = paste0(wk.dir, "/out_kernelRegVol/min2Mb_chr1_whole_maskMidSquare_gap50up_SIMvsCsnormCp_test_seed3412_nboot10_stepxgrid0.02_ygrid0.02_VUS.csv")
out.dir = paste0(wk.dir, "/out_heatmapVUS")
### OTHER SETTINGS #############################################################
out.id = "min2Mb_chr1_whole_maskMidSquare_gap50up_SIMvsCsnormCp_test_seed3412_nboot10_stepxgrid0.02_ygrid0.02"
# Which to append cell/tissue as id?
addCTIdTo = "ref" # "subj" | "ref"

# Matrix orientation; note that Z-score is calculated per row
row = "subj"
col = "ref"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(viridis)
library(ComplexHeatmap)
### FUNCTION ###################################################################
makeHMAP <- function(mx, mx.nme, clust.rw=TRUE, clust.cl=TRUE, plotTitle=out.name){
  hm <- ComplexHeatmap::Heatmap(matrix=mx, col=viridis(n=100), na_col="gray50",
                                cluster_columns=clust.cl, cluster_rows=clust.rw, 
                                row_dend_width=unit(50,"mm"), 
                                row_names_gp=gpar(fontsize=15),
                                column_title=plotTitle, 
                                column_title_gp=gpar(fontsize=8),
                                heatmap_legend_param=list(title=mx.nme))
  return(hm)
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(out.id, "_", paste0("row", row, "_col", col))

tbl <- read.csv(file=VUScsvPath, stringsAsFactors=FALSE, header=TRUE)
id.df <- strsplit(x=tbl$id, split="_subj_|_ref_")
id.df <- do.call("rbind.data.frame", c(id.df, stringsAsFactors=FALSE))
colnames(id.df) <- c("ct", "subj", "ref")
id.df$ct <- unlist(lapply(X=strsplit(x=id.df$ct, split="_", fixed=TRUE), 
                          FUN=function(x)x[3]))
tbl <- cbind.data.frame(id.df, tbl)
tbl[[addCTIdTo]] <- paste(tbl[[addCTIdTo]], tbl$ct, sep=".")
rm(id.df)

metric.v <- unique(tbl$metric)
for(metric in metric.v){
  
  eval(parse(text=paste0(
    'mx <- reshape2::acast(tbl[tbl$metric==metric,], formula=', row, '~', col,
    ', value.var="mean")'
  )))
  mx <- data.matrix(mx)
  
  out.name <- paste0(out.id, "_", metric, "_VUS")
  
  # Raw heatmap
  pdf(file=paste0(out.dir, "/", out.name, "_RAWhmap.pdf"), height=10, width=10)
  print( makeHMAP(mx=mx, mx.nme="Raw") )
  dev.off()

  # Z-score heatmap, calculated per row 
  mean.rw <- rowMeans(mx)
  sd.rw <- apply(X=mx, MARGIN=1, FUN=sd)
  z.mx <- (mx-mean.rw)/sd.rw
  
  pdf(file=paste0(out.dir, "/", out.name, "_Zhmap.pdf"), height=10, width=10)
  print( makeHMAP(mx=z.mx, mx.nme="Z-score") )
  dev.off()
  rm(mean.rw, sd.rw, z.mx)
  
  # Fold change heatmap, calculate per row, relative to first col value
  f.mx <- log2(mx/mx[,1])
  pdf(file=paste0(out.dir, "/", out.name, "_FChmap.pdf"), height=10, width=10)
  print( makeHMAP(mx=f.mx, mx.nme="Fold-change", clust.cl=FALSE) )
  dev.off()
  rm(f.mx)
  
  rm(out.name, mx); gc()
  print(paste0(metric, " done!"), quote=FALSE)
  
}

#rm(list=ls()); gc()