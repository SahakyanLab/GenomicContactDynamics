#BaseAgeScoreVsPersist - plotting the base age scores
start_time <- Sys.time()

################################################################################
#Directories

#functions to source
#lib <- "/Users/ltamon/ProjectBoard/lib"
lib <- "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for task
#objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/Plots"
objective.dir <- "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/Plots"
#BaseAgeVsPersist files
#baseAge.persist.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/out_BaseAgeVsPersist"
baseAge.persist.dir <- "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/out_BaseAgeVsPersist"
output.dir <- paste0(objective.dir, "/out_score") 

################################################################################
#Packages

library(ggplot2)
#library(foreach)

################################################################################
#Set values

#2(2Mb) or "05"(0-5 Mb), minimum gap acceptedbetween contacting bins 
#gc.v <- c(2, "05") 
gc <- 2
#chromosome
#chr.v <- c(1:22, "X")
chr <- 22
#Column names in the BaseAgeVsPersist files
Score <- "score"
Ntis <- "ntis"
plotOnly <- FALSE

################################################################################   
#Functions

source(paste0(lib, "/ggLayers_persist.R"))

myplot <- function(dta) {
  ggplot(data=as.data.frame(dta), aes(x=dta[,Ntis], y=dta[,Score]))+
    geom_violin( aes(factor(dta[,Ntis]), fill=factor(dta[,Ntis])), scale = "count" )+
    guides(fill=FALSE)+
    geom_boxplot( aes(factor(dta[,Ntis])), width=0.15 )+
    ggtitle(paste0("chr", chr, " _min", gc, "Mb"))+
    scale_y_continuous(name="Base age Score")+
    scale_x_discrete(name="COUNT (non-0 Tissues/Cell lines)", limits=c(1:21))+
    scale_fill_manual( values=c( PersistScoreColour( sort(as.numeric(unique(dta[,Ntis]))) ) ) )+
    theme.persist  
}

################################################################################  
#foreach(gc=gc.v, .inorder=TRUE) %do% {
#  foreach(chr=chr.v, .inorder=TRUE) %do% {
    if(plotOnly==FALSE){
      
      load(paste0(baseAge.persist.dir, "/chr", chr, "_min", gc, 
                  "Mb_BaseAgeVsPersist.RData"))
      df <- BASEAGE.PERSIST$POOL.AND.MX
      if(is.null(df)){
        stop(paste0("No data for Chromosome ", chr, "."))
      }
      
      scoresPerNtis.lst <-  by( df[,Score], df[,Ntis], 
                                function(x){
                                  as.vector(x)
                                } )
      scoresPerNtis.lst <- lapply( as.list(names(scoresPerNtis.lst)), 
                                   function(x){
                                     perNtis <- scoresPerNtis.lst[[x]]
                                     cllpseString <- paste(perNtis, collapse=";")
                                     spltToNum <- strsplit(cllpseString, split=";")[[1]]
                                     cbind(rep(x, length(spltToNum)), spltToNum)
                                   } )
      scoresPerNtis.df <- do.call(rbind, scoresPerNtis.lst)
      scoresPerNtis.df <- cbind(as.numeric(scoresPerNtis.df[,1]), 
                                as.numeric(scoresPerNtis.df[,2]))
      colnames(scoresPerNtis.df) <- c(Ntis, Score)
      #scoresPerNtis.df <- scoresPerNtis.df[order(scoresPerNtis.df[,Ntis]),]
      
      p <- myplot(scoresPerNtis.df)
      save(BASEAGEscore.PERSIST=scoresPerNtis.df, 
           file=paste0(output.dir, "/chr", chr, "_min", gc, 
                       "Mb_BaseAgeScoreVsPersist.RData"))
    } else {
      load(file=paste0(baseAge.persist.dir, "/chr", chr, "_min", gc, 
                       "Mb_BaseAgeVsPersist.RData"))
      p <- myplot(scoresPerNtis.df)
    }
    pdf(file=paste0(output.dir, "/chr", chr, "_min", gc, 
                       "Mb_BaseAgeScoreVsPersist.pdf"),
        height=8, width=10)
    print(p)
    dev.off()
#  }
#}
end_time <- Sys.time()
end_time-start_time   

#rm(list=ls())

