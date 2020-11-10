################################################################################
# Plot distribution of HiC signal count per tissue per chromosome
################################################################################ 
#functions to source
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
#HiC_Human21 persist files
#persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
persist.dir = "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
#output.dir = paste0(objective.dir, "/out_AnnotationVsPersist_test")   
output.dir = paste0(objective.dir, "/out_plot_HiCstrength")    

library(foreach)
library(data.table) 
library(ggplot2) 

#2(2MB gap) or "05"(0.5 MB gap), refers to minimum gap between contacting bins, 
#two points should be far enough to filter for contacts within a TAD
gcb.v = c("2", "05")

chr.v = c(1:22, "X") 

#celltis.v = c("FC")
celltis.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
              "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")

#Co = Prefrontal Cortex
#Hi = Hippocampus
#Lu = Lung (LG)
#LV = Left ventricle
#RV = Right ventricle
#Ao = Aorta
#PM = Psoas Muscle (PO)
#Pa = Pancreas
#Sp = Spleen
#Li = Liver
#SB = Small bowel
#AG = Adrenal gland
#Ov = Ovary
#B1 = Bladder
#MesC = Mesendoderm
#MSC = Mesenchymal stem cell
#NPC = Neural Progenitor cell
#TLC = Trophoblasts-like cell
#ESC = Embryonic stem cell
#FC = Fibroblast cell
#LC = Lymphoblast cell

plotOnly <- FALSE

################################################################################
source(paste0(lib, "/ggLayers_persist.R"))
source(paste0(lib, "/multiplot.R"))

myplot <- function(df=df){
  ggplot(data=as.data.frame(df), aes(x=factor(breakStartZero), y=FractionChangeHighVsLow, 
                                     group=chr, colour=factor(chr))) +
    geom_point() +
    geom_line() +
    labs(colour="Chr") +
    scale_colour_discrete(labels=as.character(chr.v)) +
    ggtitle(paste0("min", gcb, "Mb_", celltis)) +
    #xlab("HiC strength break") +
    scale_x_discrete(name="HiC strength break", breaks=breaks.plot) +
    ylab(paste0("Fraction Change in Count of High Vs Low")) +
    theme.persist #from sourced "ggLayers_persist.R"
}

gcb="2"
chr="1"
celltis="FC"

################################################################################
for(gcb in gcb.v){
  myplots.list <- list()
  for(celltis in celltis.v){
    
    if(plotOnly==FALSE){
      
      LOWHIGH.COUNT <- foreach(chr=chr.v, .inorder=TRUE, .combine=rbind) %do% {
        load(file=paste0(persist.dir, "/chr", chr, "_Persist_min", gcb, "Mb.RData"))
        N <- nrow(PERSIST.MX$hits)
        strength <- PERSIST.MX$hits[,celltis]
        rm( list=c("PERSIST.MX") )
        gc()
        #breaks from hist() covers the max value
        #h <- hist(x=strength, plot=FALSE)
        #breaks <- unique(ceiling(h$breaks))
        #breaks <- c(1, breaks[-1])
        #breaks <- breaks[-( which(breaks%in%c(0,1)) )]
        #breakend <- tail(breaks, n=1)
        #breaks <- breaks[-length(breaks)]
        
        #cut to intervals based on breaks from hist then count elements per bin
        #(x,y] left closed 
        if(chr=="X"){chr <- 23L}
        lst <- lapply(as.list(breaks), FUN=function(x){
          LH.count <- table( cut(x=strength, 
                                 breaks=c(0, x, breakend),
                                 labels=c("Low", "High"), 
                                 right=TRUE, #interval closed on the right, open on left (x,y])
                                 include.lowest=FALSE) )
          c( as.numeric(chr), x, LH.count[c("Low", "High")] )
        })
        LOWHIGH.COUNT <- do.call(rbind, lst)
        #percent change of High count relative to Low count
        add.col <- ( (LOWHIGH.COUNT[,("High")]-LOWHIGH.COUNT[,("Low")])/LOWHIGH.COUNT[,("Low")] )
        #add.col <- LOWHIGH.COUNT[,("High")]/LOWHIGH.COUNT[,("Low")]
        LOWHIGH.COUNT <- cbind(LOWHIGH.COUNT, add.col)
      }
      colnames(LOWHIGH.COUNT) <- c("chr", "breakStartZero", "Low", "High", "FractionChangeHighVsLow")
      save(LOWHIGH.COUNT, file=paste0(output.dir, "/plot_min", gcb, "Mb_", celltis, 
                                      "_chrALL_HiCstrengthBreaks.RData"))
    } else {
      load(file=paste0(output.dir, "/plot_min", gcb, "Mb_", celltis, 
                       "_chrALL_HiCstrengthBreaks.RData"))
    }
    breaks.plot <- sort( c(1, seq.int(from=0, to=ceiling( max(LOWHIGH.COUNT[,"breakStartZero"])/5 )*5, 
                           by=5)), decreasing=FALSE )
    p <- myplot(df=LOWHIGH.COUNT)
    myplots.list[[celltis]] <- p
    ggsave(file=paste0(output.dir, "/plot_min", gcb, "Mb_", celltis, 
                       "_chrALL_HiCstrengthBreaks.jpeg"), plot=p)
    ggsave(file=paste0(output.dir, "/plot_min", gcb, "Mb_", celltis, 
                       "_chrALL_HiCstrengthBreaks.pdf"), plot=p,
           units="in", width=30, height=15)
  }
  #pdf(file=paste0(output.dir, "/plot_min", gcb, 
  #                "Mb_celltisALL_chrALL_HiCstrengthBreaks.pdf"),
  #    width=30, height=15)
  #multiplot(plotlist=myplots.list, cols=7)
  #dev.off()
}

#rm(list=ls())




