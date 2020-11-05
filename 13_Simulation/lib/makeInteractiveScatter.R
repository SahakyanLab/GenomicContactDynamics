################################################################################
# Function to make interactive plot
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(plotly)
### FUNCTION ###################################################################
makeInteractiveScatter <- function(
  out.dir = "/dir",
  out.name = "plot",
  df = df$mfc,
  x = gene.x[i],
  y = gene.y[i],
  # If z is NULL, make 2D scatter
  z = gene.z[i],
  #lab.v = c(x=x, y=y, z=z),
  grp = group,
  xlim.v = NULL,
  ylim.v = NULL,
  zlim.v = NULL,
  point.lab = rownames(df),
  col.v = c("red", "blue", "yellow"),
  main = "plot"
){
  
  col.v <- sapply(X=col.v, FUN=function(x) paste0("'", x, "'"))
  col.v <- paste(col.v, collapse=",")
  
  if( !is.null(z) ){
    
    eval(parse(text=paste0(
      "p <- plot_ly(data=as.data.frame(df), x=~", x, ", y=~", y, ", z=~", z, 
      ", mode='markers', symbol=~", grp, ", symbols=1, color=~", grp, 
      ", colors=c(", col.v, "), text=~", point.lab, " )"
    )))
    p <- layout(p, scene=list(xaxis=list(title=x, range=xlim.v),
                              yaxis=list(title=y, range=ylim.v),
                              zaxis=list(title=z, range=zlim.v)),
                title=main, showlegend=TRUE)
    
  } else {
    
    eval(parse(text=paste0(
      "p <- plot_ly(data=as.data.frame(df), x=~", x, ", y=~", y, 
      ", mode='markers', symbol=~", grp, ", symbols=1, color=~", grp, 
      ", colors=c(", col.v, "), text=~", point.lab, " )"
    )))
    p <- layout(p, scene=list(xaxis=list(title=x, range=xlim.v),
                              yaxis=list(title=y, range=ylim.v)),
                title=main, showlegend=TRUE)
  }
  
  p <- add_markers(p, marker=list(size=5, opacity=0.5))
  
  htmlwidgets::saveWidget(widget=as_widget(p), 
                          file=paste0(out.dir, "/", out.name, "_plot.html"))
}
################################################################################
