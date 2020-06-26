################################################################################
# Gplot for MCL result
### FUNCTION ###################################################################
makeClustgplot <- function(MCLeq.mx=MCLOUT$equilibrium.state, coord.mx=POS.MX, 
                           title="mcl_gplot", label=node.v, v.col=NULL, v.cex=1 
){
  v.b.col <- v.col
  return(
    gplot(MCLeq.mx, gmode="graph", coord=coord.mx, main=title, cex.main=0.5,
         jitter=TRUE, displaylabels=TRUE, label=label, label.cex=0.1,
         vertex.cex=v.cex, vertex.col=v.col, vertex.border=v.b.col, edge.col="white")
  )
}
