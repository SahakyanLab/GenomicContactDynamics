############################################################################### 
#calculating fraction of human-, primate-, mammalian/ancient-specific bases
start_time <- Sys.time()

###############################################################################  
lib = "/Users/ltamon/ProjectBoard/lib"
#lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for task
objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
hpm.ntis.dir = paste0(objective.dir, "/out_ancestor")
#dir.create(paste0(objective.dir, "/out_HPMPercent"))
output.dir = paste0(objective.dir, "/out_HPMPercent") 

library(ggplot2)
library(foreach)
#install.packages("matrixStats")
library(matrixStats) #for colMedians
#install.packages("ggpubr")
library(ggpubr)

#2(2Mb) or "05"(0.5 Mb), minimum gap acceptedbetween contacting bins 
gcb = "05" #c("2", "05") 

chr.v = c(1:22, "X")

#what to plot
stat.plot.v = "means" #c("means", "medians")

#minimum allowed number of age data per bin
minAgeCountPerBin = 3000L 

plotOnly <- FALSE

#gcb="2"
#chr="21"

################################################################################ 
source(paste0(lib, "/ggLayers_persist.R"))
#source(paste0(lib, "/multiplot.R"))

myplot <- function(df=df, statname=statname, pointshp=shp){
  ggplot(data=as.data.frame(df), 
         aes(x=as.character(ntis), y=value, group=factor(HPM), colour=factor(HPM))) +
    geom_point(size=2.5, shape=shp,
               stroke=1) +
    geom_line(size=0.7) +
    annotate("text", x=ntis.uniq, y=text.loc, label=HPMPerc.NTIS$ntisCount, size=2) + 
    labs(colour="") + 
    scale_colour_manual(labels = c("Human", "Primate", "Mammal"), 
                        values=c("#c23b22", "#779ecb", "gray")) +
    ggtitle(paste0("min", gcb, "Mb_chr", chr, "_", minAgeCountPerBin, "_", statname)) +
    scale_x_discrete(name="Number (non-0 Tissues/Cell lines)",
                     limits=ntis.uniq) +
    ylab(paste0("Y")) +
    theme.persist #from sourced "ggLayers_persist.R"
}

################################################################################  
################################################################################  
myplots.list <- list()
stat.v = c("means", "medians")
for(chr in chr.v){
  
  if(plotOnly==FALSE){
    load(file=paste0(hpm.ntis.dir, "/plot_min", gcb, "Mb_chr", chr, "_", 
                     minAgeCountPerBin, "_HPMfraction_AncestorVsPersist.RData"))
    
    meanALL <- colMeans(HPM.NTIS.MX[,3:5], na.rm=TRUE)
    medianALL <- colMedians(HPM.NTIS.MX[,3:5], na.rm=TRUE)
    incl.ij.ind <- which(!is.na(HPM.NTIS.MX[,"TotalNumberOfBase"]))
    
    #count data points per ntis (for the plot)
    #sorted: increasing ntis  
    ntis.count <- table(HPM.NTIS.MX[incl.ij.ind, "ntis"])
  
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
    ntis.uniq <- as.numeric( names(ntis.count) ) 
    
    #for plot
    HPMPerc.NTIS[["ntisCount"]] <- ntis.count
    
    save(HPMPerc.NTIS, 
         file=paste0(output.dir, "/plot_min", gcb, "Mb_chr", chr, "_", 
                     minAgeCountPerBin, "_HPMPercent_AncestorVsPersist.RData"))
  } else {
    load(file=paste0(output.dir, "/plot_min", gcb, "Mb_chr", chr, "_", 
                     minAgeCountPerBin, "_HPMPercent_AncestorVsPersist.RData"))
    ntis.uniq <- sort( as.numeric(unique(HPMPerc.NTIS[[1]][,"ntis"])), decreasing=FALSE)
  }
  
  for(stat in stat.plot.v){
    if(stat=="means"){
      lab <- "mean" 
      shp <- 16
    } else {
      lab <- "median" 
      shp <- 4
    }
    #where to put ntisCount text in plot
    text.loc <- min(HPMPerc.NTIS[[stat]]["value"])-1 
    p <- myplot(df=HPMPerc.NTIS[[stat]], statname=lab, pointshp=shp)
    myplots.list[[paste0(chr, stat)]] <- p
    ggsave(file=paste0(output.dir, "/plot_min", gcb, "Mb_chr", chr, 
                       "_", minAgeCountPerBin, "_", stat,
                       "_HPMPercent_AncestorVsPersist.jpeg"),
           plot=p)
  }
}

#pdf(file=paste0(output.dir, "/plot_min", gcb, "Mb_chr", chr, 
#                "_", minAgeCountPerBin, "_", stat,
#                "_HPMPercent_AncestorVsPersist.pdf"),
#    width=13, height=8.5)
#multiplot(plotlist=myplots.list, cols=2)
#dev.off()
myplots.arranged <- ggpubr::ggarrange(plotlist=myplots.list, nrow=1, ncol=4,
                                      common.legend=TRUE, legend="bottom")
ggpubr::ggexport(myplots.arranged, width=30, height=10,
                 filename=paste0(output.dir, "/plot_min", gcb, "Mb_", minAgeCountPerBin, 
                                 "_", paste(stat.plot.v, collapse=""),
                                 "_HPMPercent_AncestorVsPersist.pdf"))

#rm(list=ls())


  