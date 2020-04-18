################################################################################
# Make Hi-Cs and Hi-Cp map of two regions
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(ggplot2)
# library(RColorBrewer)
# source(paste0(lib, "/GG_bgr.R"))
# source(paste0(lib, "/multiplot.R"))
################################################################################
hmplot <- function(df = MX, 
                   #X = bins.x,
                   #Y = bins.y,
                   out.name = paste0(chr, "_", gcb, "_", ct, "_", region),
                   scalebr.v = c(xmin=1, xmax=50, ymin=1, ymax=30),
                   outtype = "jpeg" # "png"
                   ){

  fill.lab <- list(Cs=expression(bold( "C"["s"] )),
                   Cp=expression(bold( "C"["p"] )),
                   CII=expression(bold( "C"["||"] ))
  )
  
  p.lst <- list()
  metric.v <- colnames(df)[!colnames(df)%in%c("i", "j", "h")]
  
  for(metric in metric.v){
    
    if(metric=="Cs"){
      # Blues from brewer.pal(name="Blues")
      coul <- c(brewer.pal(name="Blues", n=9)[c(2:3,5,7,9)],
                # Yellow, Red, White, Black
                "#FDC776","#9E0142", "#ffffff", "#000000")
      #brewer.pal(name="Reds", n=9)[c(5,8)]
      coul <- coul[c("1", "2", "3", "4", "5", "(5,10]", ">10", "0", "-2")%in%unique(df$Cs)]
    } else if(metric=="Cp"){
      coul <- c( colorRampPalette( rev( brewer.pal(11, "Spectral") ) )(21),
                 "#000000" )
      # +1 in case min is 0
      coul <- coul[as.character(c(1:21, -2))%in%unique(df$Cp)]
    } else if(metric=="CII"){
      coul <- c("#4292C6", "#FDC776","#9E0142", "#000000")
      coul <- coul[c("-1", "0", "1", "-2")%in%unique(df$CII)]
    } else {
      stop("Invalid type.")
    }
   
    p.lst[[metric]] <- ggplot(data=df, aes_string(x="i", y="j", z=metric)) + 
      
      geom_tile(aes_string(fill=metric, height="h")) + 
      scale_fill_manual(values=coul, na.translate=TRUE, na.value="#ffffff") + 
      labs(x=NULL, y=NULL, title=out.name,
           fill=fill.lab[[metric]]
           ) +
      bgr2 + 
      guides(fill=guide_legend(nrow=1, label.position="bottom",
                               title.position="top" #, title.hjust=0.5
      )
      ) + 
      theme(#axis.text.y=element_blank(),
        axis.text.x=element_text(face="bold", size=10, 
                                 angle=360, colour="black"),
        axis.text.y=element_text(face="bold", size=10, 
                                 angle=360, colour="black"),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=10, face="bold"),
        legend.title=element_text(size=10, face="bold"),
        legend.position="bottom",
        legend.direction="horizontal", 
        legend.box="horizontal"
        #legend.spacing.x=unit(0, 'cm'),
        #aspect.ratio=1 
      ) 
    
    if( !is.null(scalebr.v) & length(scalebr.v)==4 ){
      p.lst[[metric]] <- p.lst[[metric]] +
        geom_rect(colour="black", 
                  aes(xmin=scalebr.v["xmin"], xmax=scalebr.v["xmax"],
                      ymin=scalebr.v["ymin"], ymax=scalebr.v["ymax"]))
    }
    
  }
 
  if(outtype=="pdf"){
    pdf(file=paste0(out.dir, "/", out.name, "_sqHm.", outtype), 
        width=20, height=10)
  } else {
    jpeg(file=paste0(out.dir, "/", out.name, "_sqHm.", outtype), 
         width=20, height=10, units="in", res=1200)
  }
  
  multiplot(plotlist=p.lst, cols=3)
  dev.off()
  
}

