#BaseAgeColourRGBVsPersist
start_time <- Sys.time()

################################################################################
#Directories

#functions to source
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for task
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#BaseAgeVsPersist files
#baseAge.persist.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/out_BaseAgeVsPersist_test"
baseAge.persist.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/out_BaseAgeVsPersist"
#dir.create(paste0(objective.dir, "/out_colourRGB_test"))
#output.dir = paste0(objective.dir, "/out_colourRGB_test") 
#output.dir = paste0(objective.dir, "/out_colourRGB_test") 
output.dir = paste0(objective.dir, "/out_colourRGB") 

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
chr.v = "21" #c(1:22, "X")

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
  ggplot(data=as.data.frame(ColRGBPerNtisCount.mx), 
         aes(x=as.character(ntis), y=RGBcount, fill=factor(R1G2B3))) +
    geom_bar(position="fill", stat = "identity") +
    stat_summary(fun.data=give.n, geom="text", size=2) +
    labs(fill="") +
    scale_fill_manual(labels=c("human", "mammal", "primate"),
                      values=c("#c23b22", "gray", "#779ecb")) +
    ggtitle(paste0("min", gc, "Mb_chr", chr, "_", minAgeCountPerBin)) +
    scale_y_continuous(name="Fraction of Base Age Colour", labels=percent_format()) +
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
    ColRGBPerNtisString.lst <- by(BASEAGE.PERSIST$POOL.AND.MX[ind.incl, "colourRGB"], 
                                  BASEAGE.PERSIST$POOL.AND.MX[ind.incl, "ntis"],
                                  function(x){
                                    as.vector(x)
                                    #cllpseString <- paste(x, collapse=";")
                                  })
    ntis.uniq <- sort(as.numeric(unique(names(ColRGBPerNtisString.lst))))
    ColRGBPerNtisCount.lst <- lapply( as.list(names(ColRGBPerNtisString.lst)), 
                                      function(x){
                                        cllpseString <- paste(ColRGBPerNtisString.lst[[x]], collapse=";")
                                        total <- (nchar(cllpseString)-nchar(gsub(";", "", cllpseString)))+1
                                        red.count <- (nchar(cllpseString)-nchar(gsub("255,", "", cllpseString)))/4
                                        blue.count <- (nchar(cllpseString)-nchar(gsub(",255", "", cllpseString)))/4
                                        #total minus reds and blues = gray
                                        #grey.count <- total-( (nchar(cllpseString)-nchar(gsub("255", "", cllpseString)))/3 )
                                        grey.count <- total-(red.count+blue.count)
                                        #if(total!=(red.count+blue.count+grey.count)){
                                        #  stop("Colours do not add up.")
                                        #}  
                                        #cllpseString <- ColRGBPerNtisString.lst[[x]]
                                        #spltToRGBval <- strsplit(cllpseString, 
                                        #                         split=";")[[1]]   
                                        #total <- length(spltToRGBval)
                                        #red.count <- length(grep("255,", spltToRGBval))
                                        #blue.count <- length(grep(",255", spltToRGBval))
                                        #grey.count <- length(grep("255", spltToRGBval, 
                                        #                          invert=TRUE))
                                        #if(total!=(red.count+blue.count+grey.count)){
                                        #  stop("Colours do not add up.")
                                        #}
                                        rbind( as.numeric(c(chr,x,total,1, red.count)),
                                               as.numeric(c(chr,x,total,2, grey.count)),
                                               as.numeric(c(chr,x,total,3, blue.count)))
                                      } )
   
    ColRGBPerNtisCount.mx <- do.call(rbind, ColRGBPerNtisCount.lst)
    colnames(ColRGBPerNtisCount.mx) <- c("chr",
                                         "ntis", 
                                         "totalcount",
                                         "R1G2B3",
                                         "RGBcount")
    save(ColRGBPerNtisCount.mx, 
         file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_", minAgeCountPerBin, 
                     "ovlap_Count_BaseAgeColRGBVsPersist.RData"))
    p <- myplot(df=ColRGBPerNtisCount.mx)
    #myplots.list[[chr]] <- p
    ggsave(file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_", minAgeCountPerBin, 
                       "ovlap_Count_BaseAgeColRGBVsPersist.jpeg"),
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

#as.vector=Time difference of 0.3234868 sec - faster
#no as.vector = Time difference of 0.3288288 secs


