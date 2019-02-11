#BaseAgeScoreVsPersist - plotting the base age scores
start_time <- Sys.time()

################################################################################
#Directories

#functions to source
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for task
#objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#BaseAgeVsPersist files
#baseAge.persist.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/out_BaseAgeVsPersist_test"
baseAge.persist.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/out_BaseAgeVsPersist"
dir.create(paste0(objective.dir, "/out_score_new"))
#output.dir = paste0(objective.dir, "/out_score_test") 
output.dir = paste0(objective.dir, "/out_score_new")

################################################################################
#Packages

library(ggplot2)
library(foreach)

################################################################################
#Set values

#2(2Mb) or "05"(0-5 Mb), minimum gap acceptedbetween contacting bins 
gc.v = "2" #c("2", "05") 

#chromosome
chr.v = c(19, 21, 22) #c(1:22, "X")

#minimum allowed number of age data per bin
minAgeCountPerBin <- 1500L 

plotOnly = FALSE

#gc="2"
#chr="21"

################################################################################   
#Functions

source(paste0(lib, "/ggLayers_persist.R"))
source(paste0(lib, "/multiplot.R"))

#give data points per ntis
give.n <- function(x){ return(c(y=-0.5, label = length(x))) }

myplot <- function(df=df, statname=statname, pointshp=shp){
  ggplot(data=as.data.frame(df), 
         aes(x=as.character(ntis), y=value, group=stat)) +
    geom_point(size=2.5, shape=shp,
               stroke=1) +
    geom_line(size=0.7) +
    stat_summary(fun.data=give.n, geom="text", size=2) + 
    #labs(colour="") + 
    #scale_color_manual(labels = c("Mean", "Median"), 
    #                   values = c("#008a8c", "#ff7472")) +
    ggtitle(paste0("min", gc, "Mb_chr", chr, "_", minAgeCountPerBin, "_", statname)) +
    scale_x_discrete(name="Number (non-0 Tissues/Cell lines)",
                     limits=as.character(ntis.uniq)) +
    ylab("Base Age Score") +
    theme.persist #from sourced "ggLayers_persist.R"
}

################################################################################   
################################################################################  
foreach(gc=gc.v, .inorder=TRUE) %do% {
  myplots.list <- list()
  for(chr in chr.v){
    load(paste0(baseAge.persist.dir, "/min", gc, "Mb_chr", chr, 
                "_BaseAgeVsPersist.RData"))
    sepCounts <- t( sapply(BASEAGE.PERSIST$POOL.AND.MX[,"count"], 
                           function(x){as.numeric(strsplit(x, ";")[[1]])}) )
    ind.incl <- which(sepCounts[,1]>=minAgeCountPerBin & sepCounts[,2]>=minAgeCountPerBin)
    scoresPerNtis.lst <- by(BASEAGE.PERSIST$POOL.AND.MX[ind.incl,"score"], 
                            BASEAGE.PERSIST$POOL.AND.MX[ind.incl,"ntis"],
                            function(x) paste(x, collapse=";"))
    ntis.uniq <- sort(as.numeric(unique(names(scoresPerNtis.lst))))
    scoresPerNtis.lst <- lapply( as.list(names(scoresPerNtis.lst)), 
                                 function(x){
                                   perNtis.str <- scoresPerNtis.lst[[x]]
                                   spltToNum <- as.numeric( strsplit(perNtis.str, split=";")[[1]] )
                                   ave <- mean(spltToNum)
                                   med <- median(spltToNum)
                                   count <- length(spltToNum)
                                   #as.numeric(c(chr,x,count, ave, med))
                                   #melted mx for plotting ave and med in one graph    
                                   rbind( as.numeric(c(chr,x,count,1, ave)),
                                          as.numeric(c(chr,x,count,2, med)) )
                                 } )
    AveMedScorePerContact.mx <- do.call(rbind, scoresPerNtis.lst)
    colnames(AveMedScorePerContact.mx) <- c("chr", 
                                            "ntis", 
                                            "numberOfScores",
                                            "stat",
                                            "value")
    save(AveMedScorePerContact.mx, 
         file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_", 
                     minAgeCountPerBin, "ovlap_AveMedBaseAgeScoreVsPersist.RData"))
    for(x in 1:2){
      sbset <- subset(AveMedScorePerContact.mx,
                      AveMedScorePerContact.mx[,"stat"]==x)
      if(x=="1"){
        lab <- "mean" 
        shp <- 16
      } else {
        lab <- "median" 
        shp <- 4
      }
      p <- myplot(df=sbset, statname=lab, pointshp=shp)
      myplots.list[[paste0(x, "_", chr)]] <- p
      ggsave(file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_",lab, "_", 
                         minAgeCountPerBin, "ovlap_AveMedBaseAgeScoreVsPersist.jpeg"),
             plot=p)
    }
  }
  pdf(file=paste0(output.dir, "/plot_min", gc, "Mb_chr", paste(chr.v, collapse=""), "_",
                  minAgeCountPerBin, "ovlap_AveMedBaseAgeScoreVsPersist.pdf"),
      width=13, height=8.5)
  multiplot(plotlist=myplots.list, cols=3)
  dev.off()
}

end_time <- Sys.time()
end_time-start_time   

#rm(list=ls())

#if(plotOnly==FALSE){
  
#} else {
#  load(file=paste0(output.dir, "/min", gc, "Mb_chr", chr, 
#                   "_BaseAgeScoreVsPersist.RData"))
#}

