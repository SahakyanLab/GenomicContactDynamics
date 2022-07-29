################################################################################
# Make ggplot-based hexbin plot

# Based from make.hexbin.ggplot(): https://gist.github.com/wahalulu/1376861
## The following function cretes a hexbinplot from a dataset using ggplot
## for graphing but hexbin for doing the hexbinning.
## xvar and yvar are the x and y data
## bins is the number of bins to pass into the hexbin function (number of hexagons 
## across the x axis)
## cuts is the number of steps to use in a discrete color scale
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(hexbin)
# library(ggplot2)
# library(Hmisc)
# source(paste0(lib, "/GG_bgr.R"))
################################################################################
makeHexbinggplot <- function(xvar, 
                             yvar, 
                             bins=30, 
                             cuts=12,
                             #breaks.y=waiver(),
                             #limits.y=c(NA,NA),
                             xlab=expression( bold("C"["p"]) ),
                             ylab=bquote( bold("C"["s"]~.(affix)) ),
                             title=paste0("chrAll_", gcb, "_", ct),
                             # Number of colors should be >=cuts
                             col=viridis(cuts)) {
  
  hex <- hexbin(x=xvar, y=yvar, xbins=bins)
  gghex <- data.frame(hcell2xy(hex), count=hex@count,
                      xo=hex@xcm, yo=hex@ycm,
                      c=cut2(x=hex@count, g=cuts))

  p <- ggplot(data=gghex, aes(x=x, y=y, fill=c) ) +
    geom_hex(bins=bins) +
    scale_fill_manual(values=col) + 
    #scale_y_continuous(breaks=breaks.y, limits=limits.y) + 
    #scale_x_continuous(breaks=1:21, 
    #                   labels=as.vector(rbind(
    #                     seq(from=1, to=21, by=2),
    #                     rep("", 11)
    #                   ))[-22]
 
  #) +
    labs(title=title, x=xlab, y=ylab, fill="Count") +
    bgr2 +
    theme(legend.text=element_text(size=20, face="bold"),
          legend.title=element_text(size=25, face="bold")
          )
  
  return( list(gghex=gghex, hexplot=p) )
}
################################################################################
################################################################################
# Original code
#make.hexbin.ggplot <- function(xvar, yvar, bins, cuts) {
#  hex <- hexbin(xvar, yvar, bins)
#  gghex <- data.frame(hcell2xy(hex), count = hex@count,
#                      xo = hex@xcm, yo = hex@ycm,
#                      c = cut2(hex@count, g = cuts))
#  plot <- ggplot(gghex) +
#    geom_hex(aes(x = x, y = y, fill = c),
#             color = "black", stat = "identity") 
#  list(gghex = gghex, hexplot = plot)
#}

