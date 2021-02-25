############################################################################### 
# Heatmap
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
#library(ComplexHeatmap)
### FUNCTION ###################################################################
myheatmap <- function(mx=ELMTISSDYN.MX, colScheme=viridis::viridis(n=10), 
                      rep.group=rep.group, mx.nme="raw", cluster=FALSE, at.v=at.v
){
  
  if( grepl(x=rep.group, pattern="subfam", fixed=TRUE) ){
    h1 <- ComplexHeatmap::Heatmap(matrix=mx, col=colScheme, na_col="gray50", 
                                  cluster_columns=FALSE, cluster_rows=cluster,
                                  row_names_gp=gpar(fontsize=2),
                                  heatmap_legend_param=list(title=mx.nme, at=at.v))
    return(h1)
  } else if(rep.group=="fam"){
    h1 <- ComplexHeatmap::Heatmap(matrix=mx, col=colScheme, na_col="gray50",
                                  cluster_columns=FALSE, cluster_rows=cluster, 
                                  row_dend_width=unit(50,"mm"), 
                                  row_names_gp=gpar(fontsize=15),
                                  heatmap_legend_param=list(title=mx.nme, at=at.v))
    return(h1)
  } else {
    stop("Invalid rep.group argument.")
  }
  
}