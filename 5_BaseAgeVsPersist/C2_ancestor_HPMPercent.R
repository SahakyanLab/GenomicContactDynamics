#BaseAgeAncestorVsPersist
start_time <- Sys.time()

###############################################################################  
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for task
#objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
HPM.NTIS.dir = paste0(objective.dir, "/out_ancestor") 
#dir.create(paste0(objective.dir, "/out_ancestor_test"))
#output.dir = paste0(objective.dir, "/out_ancestor_test") 
dir.create(paste0(objective.dir, "/out_ancestor_HPMPercent"))
output.dir = paste0(objective.dir, "/out_ancestor_HPMPercent") 


library(ggplot2)
library(foreach)
library(matrixStats) #for colMedians

#2(2Mb) or "05"(0-5 Mb), minimum gap acceptedbetween contacting bins 
gc.v = "2" #c("2", "05") 

chr.v = "21" #c(19, 21, 22) #c(1:22, "X")

#minimum allowed number of age data per bin
minAgeCountPerBin <- 1000L 

#gc="2"
#chr="21"

################################################################################ 
source(paste0(lib, "/ggLayers_persist.R"))
source(paste0(lib, "/multiplot.R"))

#give data points per ntis
#give.n <- function(x){ return(c(y=-0.5, label = length(x))) }

myplot <- function(df=df, statname=statname, pointshp=shp){
  ggplot(data=as.data.frame(df), 
         aes(x=as.character(ntis), y=value, group=factor(HPM), colour=factor(HPM))) +
    geom_point(size=2.5, shape=shp,
               stroke=1) +
    geom_line(size=0.7) +
    #stat_summary(fun.data=give.n, geom="text", size=2) + 
    labs(colour="") + 
    scale_colour_manual(labels = c("Human", "Primate", "Mammal"), 
                        values=c("#c23b22", "#779ecb", "gray")) +
    ggtitle(paste0("min", gc, "Mb_chr", chr, "_", minAgeCountPerBin, "_", statname)) +
    scale_x_discrete(name="Number (non-0 Tissues/Cell lines)",
                     limits=ntis.uniq) +
    ylab("Mean/Median Fraction of Base (Percent Deviation from Overall ", lab, ")") +
    theme.persist #from sourced "ggLayers_persist.R"
}

################################################################################  
################################################################################  
foreach(gc=gc.v, .inorder=TRUE) %do% {
  for(chr in chr.v){
    load(file=paste0(HPM.NTIS.dir, "/plot_min", gc, "Mb_chr", chr, "_", 
                     minAgeCountPerBin, "_HPMfraction_AncestorVsPersist.RData"))
    
    meanALL <- colMeans(HPM.NTIS.MX[,3:5], na.rm=TRUE)
    medianALL <- colMedians(HPM.NTIS.MX[,3:5], na.rm=TRUE)
    incl.ij.ind <- which(!is.na(HPM.NTIS.MX[,"TotalNumberOfBase"]))
    
    stat.v <- c("means", "medians")
    HPM.n <- c("HumanSpecificBase", "PrimateSpecificBase", "MammalSpecificBase") 
    
    HPMPerc.NTIS <- list()
    for(stat in stat.v){
      percent.ntis <- lapply(X=as.list(seq(1:3)), FUN=function(y){
        value <- by(data=HPM.NTIS.MX[,HPM.n[y]], INDICES=HPM.NTIS.MX[,"ntis"],
                    FUN=function(x){
                      if(stat=="medians"){
                        v <- as.numeric(medianALL[y])
                        ((median(x, na.rm=TRUE)-v)/v)*100
                      } else {
                        v <- as.numeric(meanALL[y])
                        ((mean(x, na.rm=TRUE)-v)/v)*100
                      }
                    })
        data.frame( HPM=rep(y), value=as.numeric(value), ntis=names(value) )
      })
      HPMPerc.NTIS[[stat]] <- do.call(rbind, percent.ntis)
    }
    
    ntis.uniq <- sort( as.numeric(unique(HPMPerc.NTIS[[1]][,"ntis"])), decreasing=FALSE)
    
    save(HPMPerc.NTIS, 
         file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_", 
                     minAgeCountPerBin, "_HPMPercent_AncestorVsPersist.RData"))
    
    myplots.list <- list()
    for(stat in stat.v){
      if(stat=="means"){
        lab <- "mean" 
        shp <- 16
      } else {
        lab <- "median" 
        shp <- 4
      }
      p <- myplot(df=HPMPerc.NTIS[[stat]], statname=lab, pointshp=shp)
      myplots.list[[stat]] <- p
      ggsave(file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, 
                         "_", minAgeCountPerBin, "_", stat,
                         "_HPMPercent_AncestorVsPersist.jpeg"),
             plot=p)
    }
    
    pdf(file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, 
                    "_", minAgeCountPerBin, "_", stat,
                    "_HPMPercent_AncestorVsPersist.pdf"),
        width=13, height=8.5)
    multiplot(plotlist=myplots.list, cols=2)
    dev.off()
  
  }
}

#rm(list=ls())


  