################################################################################
# Wrapper function to make venn using R package venn (up to 7 sets)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(venn)
### FUNCTION ###################################################################
doVenn <- function(vennlist = vennlist,
                   # Vector of names same order as your list.
                   # If NULL, it uses names of the list.
                   labels = NULL,
                   makeVenn = TRUE,
                   saveVenndata = TRUE,
                   # If venncol="style", predefined colours will be used.
                   # Specify vector of colours to customise. 
                   venncol = "style", 
                   filename = "/out.dir/myvenn"){
  
  if( !is.null(labels) ){
    names(vennlist) <- labels
  }

  if(makeVenn==TRUE){
    pdf(file=paste0(filename, "_venn.pdf"))
    v <- venn::venn( x=vennlist, zcolor=venncol )
    dev.off()
  } else {
    v <- venn::venn( x=vennlist, zcolor=venncol )
  }
  
  #Ref: https://www.r-bloggers.com/working-with-venn-diagrams/
  # You can inspect the contents of this object with the str() function
  #str(items)
  
  # By inspecting the structure of the a object created, 
  # you notice two attributes: 1) dimnames 2) intersections
  # We can store the intersections in a new object named inters
  VENN.DATA <- attr(v, "intersections")
  
  if(saveVenndata==TRUE){
    save(VENN.DATA, file=paste0(filename, "_venn.RData"))
  }
  
}


