#Histogram to show distribution of the number of overlapping recombination rates per bin

################################################################################
#Directories

#functions to source
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#RecombinationRatesVsPersist files
#baseAge.persist.dir = paste0(objective.dir, "/out_BaseAgeVsPersist_test")
baseAge.persist.dir = paste0(objective.dir, "/out_BaseAgeVsPersist")
#dir.create(paste0(objective.dir, "/out_plot_AgeCountPerNtis_density_test"))
#output.dir = paste0(objective.dir, "/out_plot_AgeCountPerNtis_density_test")  
output.dir = paste0(objective.dir, "/out_plot_AgeCountPerNtis_density")  

################################################################################
#Packages

##CRAN
library(ggplot2)
library(foreach)

################################################################################
#Set values 

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify 
#a contact, two points should be far enough to filter for contacts within a TAD
gc.v = "2" #c("2", "05")

#write as a vector of chromosome number or X 
#no recombination data for Y and MT from Myers2005/Munsch2014
chr.v = "19" #c(1:22, "X")

plotOnly = FALSE 

################################################################################
#Functions

source(paste0(lib, "/ggLayers_persist.R"))

densplot <- function(df=df){
  ggplot(data=as.data.frame(AGECOUNT.NTIS.MX), aes(count)) +
    geom_density(position = "stack", aes(y=..count.., fill=factor(ntis))) +
    geom_vline(aes(xintercept=mean(count)),
               color="black", linetype="dashed", size=0.5) +
    labs(fill="Persistence score") +
    scale_fill_manual(values=PersistScoreColour( as.numeric(sort(unique(AGECOUNT.NTIS.MX[,"ntis"]))) ) ) +
    ggtitle(paste0("chr", chr, " _min", gc, "Mb")) +
    xlab("Number of bases with age data per bin") +
    #include number of datapoints per ntis
    ylab(paste0("Density (count)")) +
    theme.persist #from sourced "ggLayers_persist.R"
}
#
################################################################################
foreach(chr=chr.v, .inorder=TRUE) %do% {
  for(gc in gc.v){
    
    if(plotOnly==FALSE){
      
      load(paste0(baseAge.persist.dir, "/min", gc, "Mb_chr", chr, 
                  "_BaseAgeVsPersist.RData"))
      
      ntis.uniq <- sort(unique(as.numeric(BASEAGE.PERSIST$POOL.AND.MX[,"ntis"])), 
                        decreasing=FALSE)
      
      bins.uniq.ntis <- by(BASEAGE.PERSIST$POOL.AND.MX[, c("i", "j")],
         BASEAGE.PERSIST$POOL.AND.MX[,"ntis"], function(x){
           unique(unlist(x))
         })
      
      lst <- lapply(as.list(names(bins.uniq.ntis)), function(x){
        perNtis <- bins.uniq.ntis[[x]]
        mtched <- match(perNtis, BASEAGE.PERSIST$AGE.BIN.MX[,c("bin")])
        count.ntis <- BASEAGE.PERSIST$AGE.BIN.MX[mtched,"count"]
        cbind(ntis=rep(as.numeric(x)), 
              bin=as.numeric(perNtis), 
              count=as.numeric(count.ntis))
      })
    
      AGECOUNT.NTIS.MX <- do.call(rbind, lst)
      save(AGECOUNT.NTIS.MX, file=paste0(output.dir,"/plot_min", gc, "Mb_chr", 
                                         chr, "_AgeCountPerNtis.RData"))
    } else {
      load(file=paste0(output.dir,"/plot_min", gc, "Mb_chr", 
                       chr, "_AgeCountPerNtis.RData"))
    }
    
    #N <- nrow(AGECOUNT.NTIS.MX)
    densplot(df=AGECOUNT.NTIS.MX)
    ggsave(file=paste0(output.dir, "/dens_each_min", gc, "Mb_chr", chr, 
                       "_AgeCountPerBin.jpeg"))
  }
}

#rm(list=ls())

