#Generate plots for recombination rates~persistent contacts 

################################################################################ 
#Set directories

#functions to source
#lib <- "/Users/ltamon/ProjectBoard/lib"
lib <- "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/2_RecombinationRatesVsPersist"
objective.dir <- "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/2_RecombinationRatesVsPersist"
#RecombinationRatesVsPersist files
recom.persist.dir <- paste0(objective.dir, "/out_RecomVsPersist")
#dir.create(paste0(objective.dir, "/out_makeBoxplot"))
output.dir <- paste0(objective.dir, "/out_makeBoxplot")  

################################################################################  
#Packages
##CRAN
library(ggplot2)
library(foreach)

################################################################################  
#Set values 

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify 
#a contact, two points should be far enough to filter for contacts within a TAD
gc.v <- c(2, "05")

#write as a vector of chromosome number or X 
#no recombination data for Y and MT from Myers2005/Munsch2014
chr.v <- c(1:22, "X")

plotOnly <- FALSE 

#sel.plot.v <- c("MEANrb", "MEDIANrb", "POOL")
#sel.stat.v <- c("MEANijP", "MEDIANijP", "SDEVijP", "MINijP", "MAXijP",
#                "MEANij", "MEDIANij", "SDEVij", "MINij", "MAXij", "DIFFij")

################################################################################   
#Functions

#source(paste0(lib, "/UTL_NCPU.R"))
source(paste0(lib, "/ggLayers_persist.R"))

#give data points per box in boxplot
give.n <- function(x){ return(c(y=-0.5, label = length(x))) }

myplot <- function(df=df){
  ggplot(data=as.data.frame(df), aes(x=ntis, y=df[,2])) +
    geom_boxplot(aes(fill=factor(ntis))) +
    guides(fill=FALSE) + 
    stat_summary(fun.data=give.n, geom="text", size=2) + 
    ggtitle(paste0("chr", chr, " _min", gc, "Mb", "_", plot, "_", stat)) +
    scale_x_discrete(name="COUNT (non-0 Tissues/Cell lines)", 
                     limits=c(1:21), drop=FALSE) +
    scale_fill_manual(values=PersistScoreColour( as.numeric(sort(unique(df[,"ntis"]))) ) ) +
    scale_y_continuous(name=paste0(stat, " RecombinationRate(cM/Mb)")) +
    theme.persist #from sourced "ggLayers_persist.R"
}

################################################################################ 
################################################################################   
#nested foreach to get data for each chromosome for each gap 
foreach(gc=gc.v, .inorder=TRUE) %do% {
  foreach(chr=chr.v, .inorder=TRUE) %do% { 
    
    if(plotOnly==FALSE){
      load(paste0(recom.persist.dir, "/chr", chr, "_min", gc, "Mb_RecomVsPersist.RData"))
      
      ntis   <- RECOM.PERSIST$NTIS
      toPlot <- RECOM.PERSIST$BINVAL
      toPlot[["POOL"]] <- RECOM.PERSIST$POOL
      
      plot.v <- c("MEANrb", "MEDIANrb", "POOL")
      
      foreach(plot=plot.v, .inorder=TRUE) %do% {
        
        if (plot=="POOL"){
          stat.v <- c("MEANijP", "MEDIANijP", "SDEVijP", "MINijP", "MAXijP")
        } else {
          stat.v <- c("MEANij", "MEDIANij", "SDEVij", "MINij", "MAXij", "DIFFij")
        }
        
        dropped.rows.NA <- c()
  
        foreach(stat=stat.v, .inorder=TRUE) %do% {
          
          #subset data to be plotted
          df <- as.data.frame( cbind(ntis, as.numeric(toPlot[[plot]][,stat])) )
          colnames(df) <- c("ntis", paste0(stat))
          
          #remove rows with NAs
          NA.rows      <- apply(df, 1, function(x){any(is.na(x))})
          NA.rows.ind  <- which(NA.rows==TRUE)
          drop.rows.NA <- sum(NA.rows) #number of dropped rows (with NA)
          if(drop.rows.NA!=0){
            df <- df[-c(NA.rows.ind),]
          }
          
          dropped.rows.NA[stat] <- drop.rows.NA
          
          save(df, 
               file=paste0(output.dir,"/bp_min", gc, "Mb_", plot, "_", stat, 
                           "_chr", chr, "_RecomVsPersist.Rdata"))
          
          p <- myplot(df=df)
          ggsave(file=paste0(output.dir, "/bp_min", gc, "Mb_", plot, "_", stat, 
                             "_chr", chr, "_RecomVsPersist.png"),
                 p)
          
          #pdf(file=paste0(output.dir, "/bp_min", gc, "Mb_", plot, "_", stat, 
          #                "_chr", chr, "_RecomVsPersist.pdf"), height=8, width=8)
          #print(p)
          #dev.off()
          
        } #stat
        
        if(sum(dropped.rows.NA)!=0){
          write.table(as.data.frame(dropped.rows.NA), 
                file=paste0(output.dir,"/bp_min", gc, "Mb_", 
                            plot, "_chr", chr, "_RecomVsPersist_SumDroppedNAs.txt"),
                col.names=FALSE, quote=FALSE)
        }
        
      } #plot
      
    } else {
      foreach(plot=sel.plot.v, .inorder=TRUE) %do% {
        foreach(stat=sel.stat.v, .inorder=TRUE) %do% {
          load(file=paste0(output.dir,"/bp_min", gc, "Mb_", plot, "_", stat, 
                           "_chr", chr, "_RecomVsPersist.Rdata"))
          
          p <- myplot(df=df)
          ggsave(file=paste0(output.dir, "/bp_min", gc, "Mb_", plot, "_", stat, 
                             "_chr", chr, "_RecomVsPersist.png"),
                 p)
          
          #pdf(file=paste0(output.dir, "/bp_min", gc, "Mb_", plot, "_", stat, 
          #                "_chr", chr, "_RecomVsPersist.pdf"), height=8, width=8)
          #print(p)
          #dev.off()
          
        } # sel.stat.v
      } #sel.plot.v
    } #else
  } #end bracket for chromosome foreach
} #end bracket for gap foreach  

#rm(list=ls())





