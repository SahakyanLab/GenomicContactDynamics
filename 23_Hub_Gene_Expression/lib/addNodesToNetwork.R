################################################################################
# Add nodes to network e.g., nodes with feature 
### FUNCTION ###################################################################
addNodesToNetwork <- function(mark.df, NETWRK, edgelen.mult){
  
  add <- unique( c(mark.df[,"bin.from"],mark.df[,"bin.to"]) )
  base <- c(1,NETWRK$edges$to[is.na(NETWRK$edges$label)])
  add <- setdiff(add, base)
  ubins <- sort(unique(c(add, base)), decreasing=FALSE)
  dot.sample <- NETWRK$nodes[NETWRK$nodes$shape=="dot",][1,]
  edge.sample <- NETWRK$edges[is.na(NETWRK$edges$label),][1,]
  NETWRK$nodes <- rbind(NETWRK$nodes,
                        data.frame(id=add, label=NA, title=paste0("Bin: ", add), 
                                   shape="square", 
                                   size=dot.sample$size, color=edge.sample$color, 
                                   borderWidth=dot.sample$borderWidth, font.size=dot.sample$font.size, #font=dot.sample$font, 
                                   stringsAsFactors=FALSE)
  )
  rm(dot.sample)
  
  NETWRK$nodes <- NETWRK$nodes[order(NETWRK$nodes$id),]
  rm(base); gc()
  
  #-------------------
  
  # Contact edges should remain the same
  # Combine contact edges to the updated edges for non-contacts (with new nodes)
  nonij.sample <- NETWRK$edges[is.na(NETWRK$edges$label),][1,]
  NETWRK$edges <- rbind(NETWRK$edges[!is.na(NETWRK$edges$label),],
                        data.frame(#id=NA, 
                          label=NA, title=NA, from=ubins[-(length(ubins))], 
                          to=ubins[-1], length=(diff(ubins)-1)*edgelen.mult,  
                          width=nonij.sample$width, color=nonij.sample$color, 
                          arrows=nonij.sample$arrows, font.size=NA, #font=NA, 
                          stringsAsFactors=FALSE)
  )
  rm(nonij.sample)
  
  NETWRK$edges <- NETWRK$edges[order(NETWRK$edges$from),]
  print("Nodes from markfile added.", quote=FALSE)
  
  return(NETWRK)
  
}