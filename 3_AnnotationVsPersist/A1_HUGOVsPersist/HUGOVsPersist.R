#Make HUGOnamesVsPersist plot using AnnotationVsPersist files

#Take the HUGOnames per persistent score (ntis which stands for number of 
#tissues present)
#Depending on the objective, one might need to differentiate between RefSeq
#categories of the annotation (e.g. NM_ vs. NR_)
#Below are approaches with and without separation of RefSeq categories (more
#details provided below)

#To see the RefSeq categories in your annotation file, use the function
#RefSeqCat in lib/RefSeqCat.R

#REMINDERS BEFORE RUNNING
#Here, one cannot use annotation files with other RefSeq categories besides 
#NM_ and NR_ 
#Does not matter for now because hg19 and hg38 annotation files 
#only have these two

#TO DO
#further stratify NR_? - write separate script

#include more data in the saved object for plotting? now, the object just 
#gives the HUGOname then the ntis
################################################################################ 
#Set directories

#functions to source
lib <- "/Users/ltamon/ProjectBoard/lib"
#main directory for the task
objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist/A1_HUGOVsPersist"
#annotation file 
annoFile <- "/Users/ltamon/Database/ucsc_tables/hg19anno"
#AnnotationVsPersistfiles 
ANNO.PERSIST.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist/out_AnnotationVsPersist"
#dir.create(paste0(objective.dir, "/out_HUGOVsPersist"))
output.dir <- paste0(objective.dir, "/out_HUGOVsPersist")  

################################################################################  
#Packages
##CRAN
#install.packages("foreach") 
library(foreach)
#install.packages("data.table") 
library(data.table) #for fread and rbindlist
#install.packages("ggplot2")
library(ggplot2)

################################################################################ 
#Set values

NCPU <- 1

#genome version
genome.ver <- "hg19"

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify a
#contact, two points should be far enough to filter for contacts within a TAD
gc.v <- c(2, "05")

#chromosomes
chr.v <- c("ALL") #c(1:22, "X", "ALL") 

#column name for persistence score (ntis, which stands for number of tissues)
#in AnnotationVsPersist files
ntis <- "ntis"

#column name of HUGO names in AnnotationVsPersist files
name2 <- "name2"

#column name of annotation name (starts with NM or NR) 
#in annotation file and in AnnotationVsPersist files
name <- "name"

#check in AnnotationVsPersist files
#minimum overlap used when intersecting annotations to persistent contacts
olap.min <- 2500 #or 1 (default)

plotOnly <- TRUE #TRUE or FALSE

################################################################################
#Functions

source(paste0(lib, "/UTL_NCPU.R"))
source(paste0(lib, "/RefSeqCat.R")) 
source(paste0(lib, "/UCSCTableReadFilter.R"))
source(paste0(lib, "/ggLayers_persist.R"))

#General format of bar chart
myplot <- function(data){
  p <- ggplot(data=data, aes(ntis)) + 
    #plot title
    ggtitle(paste0("chr", chr, " _min", gc, "Mb_", genome.ver)) +
    scale_x_discrete(name="COUNT (non-0/total Tissues)", limits=c(1:21), 
                     drop=FALSE) +
    #scale_y_discrete(name="Number of HUGO genes (unique)", limits=c(0:6), 
    #                 drop=FALSE) +
    ylab("COUNT (HUGO names)") +
    theme.persist #from sourced "ggLayers_persist.R"
}

################################################################################ 
if(plotOnly==FALSE){

#to check for RefSeq categories in the annotation file, use the RefSeqCat function
  RefSeqCat(filePath = "/Users/ltamon/Database/ucsc_tables/hg19anno",
            genomeVersion = "hg19",
            accessionCol = "name")
  
  #BAR CHARTS - HUGOVsPersist
  foreach(gc=gc.v, .inorder=FALSE) %op% {
    foreach(chr=chr.v, .inorder=FALSE) %op% {
      
      #load ANNO.PERSIST.MX table
      load(paste0(ANNO.PERSIST.dir, "/", "chr", chr, "_min", gc, "Mb_", 
                  genome.ver, "_AnnoVsPersist_olapMin", olap.min, ".RData"))
      
      #------------------------------------------------
      #QUERY ANNO.PERSIST.MX, determine persistent score (ntis) for annotations 
      #then plot in bar chart
      #2 cases will be done:
      #1. DO NOT SEPARATE annotations into RefSeq categories, NM and NR in this case 
      #because for the functional annotation, we just need the HUGO names
      #2. SEPARATE annotations into RefSeq categories
      #because some HUGO names have both NM and NR variants and we might just 
      #focus on either of the classes
      #unique(HUGO names)
      #store HUGOnames~persist(ntis) matrices for functional annotation
      
      #separate ANNO.PERSIST.MX into classes
      ANNO.PERSIST.NM <- UCSCTableReadFilter(Table = ANNO.PERSIST.MX,
                                             Filtering.Scheme = "Human.Nuclear.mRNA")
      
      ANNO.PERSIST.NR <- UCSCTableReadFilter(Table = ANNO.PERSIST.MX,
                                             Filtering.Scheme = "Human.ncRNA")
      
      #classifies annotation in ANNO.PERSIST.MX according to ntis
      #one annotation can be assigned to >1 ntis
      lst    <- list()
      lst.NM <- list()
      lst.NR <- list()
      for(i in 1:21){
        ind         <- grep(i, ANNO.PERSIST.MX[,ntis])
        lst[[i]]    <- cbind(ANNO.PERSIST.MX[ind,], 
                             ntis.plot=rep(i, length(ind)))
        
        ind.NM      <- grep(i, ANNO.PERSIST.NM[,ntis])
        lst.NM[[i]] <- cbind(ANNO.PERSIST.NM[ind.NM,], 
                             ntis.plot=rep(i, length(ind.NM)))
        
        ind.NR      <- grep(i, ANNO.PERSIST.NR[,ntis])
        lst.NR[[i]] <- cbind(ANNO.PERSIST.NR[ind.NR,], 
                             ntis.plot=rep(i, length(ind.NR)))
      }
      
      #per ntis, take unique(HUGOnames)
      #do these for NM only, NR only and combined (no separation)
      HUGO.NTIS.list <- lapply(list(ALL=lst, NM=lst.NM, NR=lst.NR), function(y){
        lst.a  <- lapply(y, function(x) {
          lst1 <- by(x, x[,name2], function(x) unique(x))
          HUGO <- names(lst1)
          cbind(ntis=rep(x[1,"ntis.plot"], length(HUGO)), HUGO=HUGO)
        } ) 
        as.data.frame(do.call(rbind, lst.a))
      } )
      
      foreach (data=c("ALL", "NMNR"), .inorder=FALSE) %op% {
        
        pdf(file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                        "_", data, "HUGOVsPersist.pdf"), height=8, width=8)
        
        if(data=="ALL"){
          
          #save objects for the bar charts
          #no separation of annotation types (NM, NR)
          HUGO.NTIS <- HUGO.NTIS.list$ALL
          save(HUGO.NTIS, file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", 
                                      genome.ver, "_", data, "HUGOVsPersist.RData"))
          #bar chart
          p <- myplot(as.data.frame(HUGO.NTIS)) + geom_bar(stat="count")
          
        } else {
          
          #with separation of annotation types (NM, NR)
          #combine NM and NR tables differentiating between the two by 
          #adding column indicating if NM or NR
          HUGO.NTIS.NM <- cbind(Grp=rep("NM", nrow(HUGO.NTIS.list$NM)), 
                                HUGO.NTIS.list$NM)
          HUGO.NTIS.NR <- cbind(Grp=rep("NR", nrow(HUGO.NTIS.list$NR)), 
                                HUGO.NTIS.list$NR)
          HUGO.NTIS.NMNR <- rbind(HUGO.NTIS.NM, HUGO.NTIS.NR)
          
          save(HUGO.NTIS.NMNR, file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", 
                                           genome.ver, "_", data, "HUGOVsPersist.RData"))
          #bar chart
          p <- myplot(as.data.frame(HUGO.NTIS.NMNR)) + 
            geom_bar(stat="count", position="dodge2", aes(fill=Grp)) 
        }
        
        #for saving as pdf
        print(p)
        dev.off()
        
        ggsave(file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                           "_", data, "HUGOVsPersist.png"))
        
      } #end foreach data for plotting
      
    } #end foreach chr
  } #end foreach gc
  
} else { #else statement for plotOnly
  
  #load objects
  foreach(gc=gc.v, .inorder=FALSE) %op% {
    foreach(chr=chr.v, .inorder=FALSE) %op% {
      foreach (data=c("ALL", "NMNR"), .inorder=FALSE) %op% {
        
        pdf(file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                        "_", data, "HUGOVsPersist.pdf"), height=8, width=8)
        
        if(data=="ALL"){
        
          load(paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                      "_", data, "HUGOVsPersist.Rdata"))
          p <- myplot(as.data.frame(HUGO.NTIS)) + 
            geom_bar(stat="count")
          
        } else {
          
          load(paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                      "_", data, "HUGOVsPersist.Rdata"))
          p <- myplot(as.data.frame(HUGO.NTIS.NMNR)) + 
            geom_bar(stat="count", position="dodge2", aes(fill=Grp))
          
        }
        
        #for saving as pdf
        print(p)
        dev.off()
        
        ggsave(file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", 
                           genome.ver, "_", data, "HUGOVsPersist.png"))
        
      } #end foreach data for plotting
    }
  }
  
} #end else statement of plotOnly

#rm(list=ls())
