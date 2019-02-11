#BaseAgeColourRGBVsPersist
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
#dir.create(paste0(objective.dir, "/out_ancestor_test"))
#output.dir = paste0(objective.dir, "/out_ancestor_test") 
#dir.create(paste0(objective.dir, "/out_ancestor"))
output.dir = paste0(objective.dir, "/out_ancestor") 

################################################################################
#Packages

library(ggplot2)
library(foreach)
library(scales)

################################################################################
#Set values

#2(2Mb) or "05"(0-5 Mb), minimum gap acceptedbetween contacting bins 
gc.v = "2" #c("2", "05") 

#chromosome
chr.v = "21" #c(19, 21, 22) #c(1:22, "X")

#minimum allowed number of age data per bin
minAgeCountPerBin <- 1500L 

#plotOnly = FALSE

#gc="2"
#chr="21"

################################################################################   
#Functions

source(paste0(lib, "/ggLayers_persist.R"))
#source(paste0(lib, "/multiplot.R"))

#give data points per ntis
give.n <- function(x){ return(c(y=0, label = length(x))) }

myplot <- function(df=df){
  ggplot(data=as.data.frame(df), 
         aes(x=as.character(ntis), y=ancestorcount, fill=factor(ancestor))) +
    geom_bar(position="fill", stat = "identity") +
    stat_summary(fun.data=give.n, geom="text", size=2) + 
    ggtitle(paste0("min", gc, "Mb_chr", chr, "_", minAgeCountPerBin)) +
    labs(fill="") +
    scale_fill_manual(labels=c("human", "mammal", "primate"),
                      values=c("#c23b22", "gray", "#779ecb")) +
    scale_y_continuous(name="Fraction of Group", labels=percent_format()) +
    scale_x_discrete(name="Number (non-0 Tissues/Cell lines)",
                     limits=as.character(ntis.uniq)) +
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
    AncestorPerNtisString.lst <- by(BASEAGE.PERSIST$POOL.AND.MX[ind.incl, "ancestor"], 
                                  BASEAGE.PERSIST$POOL.AND.MX[ind.incl, "ntis"],
                                  function(x){
                                    as.vector(x)
                                    #cllpsestring <- paste(x, collapse=";")
                                  })
    ntis.uniq <- sort(as.numeric(unique(names(AncestorPerNtisString.lst))))
    AncestorPerNtisCount.lst <- lapply( as.list(names(AncestorPerNtisString.lst)), 
                                      function(x){
                                        cllpseString <- paste(AncestorPerNtisString.lst[[x]],
                                                              collapse=";")
                                        spltToAnc <- strsplit(cllpseString, split=";")[[1]]
                                        total <- length(spltToAnc)
                                        #counts mammals using mouse
                                        mammal.count <- length(grep("Mmus", spltToAnc, 
                                                                    ignore.case=TRUE))
                                        spltToAnc <- gsub("h|s|a|p|\\-|\\[|\\]|\\d", "", 
                                                          spltToAnc, ignore.case=TRUE)
                                        human.count <- length(which((nchar(spltToAnc))==0))
                                        primate.count <- total-(human.count+mammal.count)
                                        rbind( as.numeric(c(chr,x,total,1, human.count)),
                                               as.numeric(c(chr,x,total,2, mammal.count)),
                                               as.numeric(c(chr,x,total,3, primate.count)))
                                      } )
    AncestorPerNtisCount.mx <- do.call(rbind, AncestorPerNtisCount.lst)
    colnames(AncestorPerNtisCount.mx) <- c("chr",
                                         "ntis", 
                                         "totalcount",
                                         "ancestor",
                                         "ancestorcount")
    save(AncestorPerNtisCount.mx, 
         file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_", 
                     minAgeCountPerBin, "ovlap_Count_AncestorVsPersist.RData"))
    p <- myplot(df=AncestorPerNtisCount.mx)
    #myplots.list[[chr]] <- p
    ggsave(file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_",
                       minAgeCountPerBin, "ovlap_Count_AncestorVsPersist.jpeg"),
           plot=p)
  }
  #jpeg(file=paste0(output.dir, "/min", gc, "Mb_chr", paste(chr.v, collapse=""),
  #                 "_Count_BaseAgeColRGBVsPersist.jpeg"), width=842, height=595)
  #multiplot(plotlist=myplots.list, cols=4)
  #dev.off()
}
#rm(list=ls())
end_time <- Sys.time()
end_time-start_time   
