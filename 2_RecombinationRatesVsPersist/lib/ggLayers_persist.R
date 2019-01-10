################################################################################
PersistScoreColour <- function(PersistScores = ntisUsed){
  hexColour <- lapply(PersistScores, function(x) {
    rgb.v <- c(255,255,255)-(12*x)
    rgb(rgb.v[1], rgb.v[2], rgb.v[3], max=255)
  } )
   unlist(hexColour)
}

################################################################################
theme.persist <- theme(panel.grid.major=element_blank(),       
                       panel.grid.minor=element_blank(),
                       #modify border lines of plot
                       panel.background=element_rect(colour="black", 
                                                     size=1, fill=NA),
                       #main title
                       plot.title=element_text(face="bold", color="black", 
                                               size=14, hjust=0.5),
                       #x-axis title 
                       axis.title.x=element_text(size=13),
                       #y-axis title
                       axis.title.y=element_text(size=13),
                       #all text within theme function
                       text=element_text(color="black"),
                       #modify axis tick labels
                       axis.text.x = element_text(face="bold", size=11, angle=360),
                       axis.text.y = element_text(face="bold", size=11, angle=360)
                 )
