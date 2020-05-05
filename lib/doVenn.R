################################################################################
#make venn (up to 7 sets)
#input is a list of objects
################################################################################
library(venn)
library(gplots)
################################################################################
################################################################################
doVenn <- function(vennlist = VENN.OBJ.LIST,
                   #or a vector of names same order as your list
                   #if NULL, it uses names on the list
                   labels = NULL,
                   makeVenn = TRUE,
                   saveVenndata = TRUE,
                   venncol = "style", #or specifiy a vector
                   filename = "myvenn"){
  
  if(!is.null(labels)){
    names(vennlist) <- labels
  }
  
  if(makeVenn==TRUE){
    #make and save venn in a pdf
    pdf(file=paste0(filename, ".pdf"))
    v <- venn::venn( x=VENN.OBJ.LIST, zcolor="style" )
    dev.off()
  } else {
    v <- venn::venn( x=VENN.OBJ.LIST, zcolor="style" )
  }
  
  #save underlying venn data
  items <- gplots::venn(VENN.OBJ.LIST, show.plot=FALSE)
  
  #Ref: https://www.r-bloggers.com/working-with-venn-diagrams/
  # You can inspect the contents of this object with the str() function
  #str(items)
  
  # By inspecting the structure of the a object created, 
  # you notice two attributes: 1) dimnames 2) intersections
  # We can store the intersections in a new object named inters
  VENN.DATA <- attr(items, "intersections")
  
  if(saveVenndata==TRUE){
    save(VENN.DATA, file=paste0(filename, ".RData"))
  }
  
}


