#Make HUGOnames~NTIS plot using anno~persist
#used hg19 here

#Reminders before running
#the code can't use annotation files with other annotation types besides NM and NR (does not matter for now
#because hg19 and hg38 annotation files only have NM and NR)

#TO DO
#further stratify NR? - write separate script
#polish UCSC outsource script
#include more data in the saved object for plotting? now, the object just gives the HUGOname then the ntis
############################################################################################################ 
#Set directories

#functions to outsource
source.dir       <- "/Users/ltamon/ProjectBoard/outsource_R"
#main directory for the task
objective.dir    <- "/Users/ltamon/ProjectBoard/HiC_contacts_dim/annoGenome"
#annotation file
annoFile.path    <- "/Users/ltamon/Database/ucsc_tables"
#annotation~persistent contact files 
ANNO.PERSIST.dir <- paste0(objective.dir, "/out_anno_persist")

#Create directories 

#for outputs (i.e. ntis~HUGOname bar charts and corresponding objects)
#dir.create(paste0(objective.dir, "/out_HUGOname_ntis"))
output.dir <- paste0(objective.dir, "/out_HUGOname_ntis")  #check if it exists already

############################################################################################################  
#Packages
##CRAN
#install.packages("foreach") 
library(foreach)
#install.packages("data.table") 
library(data.table) #for fread and rbindlist
#install.packages("ggplot2")
library(ggplot2)

############################################################################################################ 
#Set values

NCPU          <- 1

#genome version
genome.ver    <- "hg19"
anno.filename <- "hg19anno" #or "hg38anno"

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gc.v          <- c(2, "05")

#no recombination data for Y and MT from Myers2005/Munsch2014
chr.v         <- c("ALL") #c(1:22, "X", "ALL") 

#column name for number of tissues(ntis) of persistent contacts in
#annotation~persistent contact files
ntis          <- "ntis"
#column name of HUGO names in in annotation~persistent contact files
name2         <- "name2"
#column name of annotation name (starts with NM or NR) 
#in annotation file and annotation~persistent contact files
name          <- "name"

#check annotation~persistent contact files
#minimum overlap used when intersecting annotations to persistent contacts
olap.min       <- 2500 #or 1 (default)

plotOnly      <- FALSE #TRUE or FALSE
############################################################################################################
#Functions

source(paste0(source.dir, "/UTL_NCPU.R"))
source(paste0(source.dir, "/UCSCTableReadFilter_edit.R"))
source(paste0(source.dir, "/ggLayers_persist.R"))

#General format of bar chart
myplot <- function(data){
  p <- ggplot(data=data, aes(x=ntis)) + 
    #plot title
    ggtitle(paste0("chr", chr, " _min", gc, "Mb_", genome.ver)) +
    scale_x_discrete(name="COUNT (non-0/total Tissues)", limits=c(1:21), drop=FALSE) +
    #scale_y_discrete(name="Number of HUGO genes (unique)", limits=c(0:6), drop=FALSE) +
    ylab("COUNT (genes)") +
    theme.persist #from sourced "ggLayers_persist.R"
}

############################################################################################################ 
if(plotOnly==FALSE){
  
  #CHECK classes in annotation file 
  #load annotation file
  anno.file <- fread(file=paste0(annoFile.path, "/", anno.filename), header=TRUE, data.table=FALSE)
  
  anno.type <- unique( sapply( as.vector(anno.file[,name]), simplify=TRUE, USE.NAMES=FALSE, 
                               function(x){
                                 b <- strsplit(x, split="_")
                                 strsplit(unlist(b), split="\t")[[1]]
                               } ) )
  
  state <- paste(anno.type, collapse=", ")
  print(paste0("Annotation types in ", genome.ver, " annotation file are ", state, "."))
  #Classes are only NM and NR. (same for hg38)
  
  #BAR CHARTS - HUGOnames~ntis
  #generate ntis~HUGOname dataframe (each for separation and no separation of NM and NR annotations) then
  #make bar charts
  foreach(gc=gc.v, .inorder=FALSE) %op% {
    foreach(chr=chr.v, .inorder=FALSE) %op% {
      
      #load ANNO.PERSIST.MX table
      load(paste0(ANNO.PERSIST.dir, "/", "chr", chr, "_min", gc, "Mb_", genome.ver, "_anno_persist_olapMin", olap.min, ".RData"))
      
      #------------------------------------------------
      #QUERY ANNO.PERSIST.MX, associate annotations with NTIS then plot in bar chart
      #2 cases:
      #1. DO NOT SEPARATE annotations into NM and NR (only two categories in annotation file name, see checking above)
      #because for the GO and pathway analyses, we just need the names
      #2. SEPARATE annotations into NM and NR
      #because some HUGO names have both NM and NR variants and we might just focus on either of the classes
      #unique(HUGO names)
      #store HUGO names for gene ontology and pathway analysis
      
      #separate ANNO.PERSIST.MX into classes
      ANNO.PERSIST.NM <- UCSCTableReadFilter( Table = ANNO.PERSIST.MX,
                                              Filtering.Scheme = "Human.Nuclear.mRNA" )
      
      ANNO.PERSIST.NR <- UCSCTableReadFilter( Table = ANNO.PERSIST.MX,
                                              Filtering.Scheme = "Human.ncRNA" )
      
      #checkpoint
      if((nrow(ANNO.PERSIST.NM)+nrow(ANNO.PERSIST.NR))!=nrow(ANNO.PERSIST.MX)){
        stop(paste0("Problem with classifying annotations into types (", state, ")."))
      }
      
      #classifies entries/anno in ANNO.PERSIST.MX into NTIS
      #one anno/entry can be assigned to >1 NTIS
      lst    <- list()
      lst.NM <- list()
      lst.NR <- list()
      for(i in 1:21){
        ind         <- grep(i, ANNO.PERSIST.MX[,ntis])
        lst[[i]]    <- cbind(ANNO.PERSIST.MX[ind,], ntis.plot=rep(i, length(ind)))
        
        ind.NM      <- grep(i, ANNO.PERSIST.NM[,ntis])
        lst.NM[[i]] <- cbind(ANNO.PERSIST.NM[ind.NM,], ntis.plot=rep(i, length(ind.NM)))
        
        ind.NR      <- grep(i, ANNO.PERSIST.NR[,ntis])
        lst.NR[[i]] <- cbind(ANNO.PERSIST.NR[ind.NR,], ntis.plot=rep(i, length(ind.NR)))
      }
      
      #names(lst)    <- c(1:21)
      #names(lst.NM) <- c(1:21)
      #names(lst.NR) <- c(1:21)
      
      #per NTIS, take unique HUGO names
      #do these for NM only, NR only and combined
      HUGO.NTIS.list <- lapply(list(ALL=lst, NM=lst.NM, NR=lst.NR), function(y){
        lst.a  <- lapply(y, function(x) {
          lst1 <- by(x, x[,name2], function(x) unique(x))
          HUGO <- names(lst1)
          cbind(ntis=rep(x[1,"ntis.plot"], length(HUGO)), hg19HUGOname=HUGO)
        } ) 
        as.data.frame(do.call(rbind, lst.a))
      } )
      
      foreach (data=c("ALL", "NMNR"), .inorder=FALSE) %op% {
        
        pdf(file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                        "_", data, "HUGO_ntis.pdf"), height=8, width=8)
        
        if(data=="ALL"){
          
          #save objects for the bar charts
          #no separation of annotation types (NM, NR)
          HUGO.NTIS <- HUGO.NTIS.list$ALL
          save(HUGO.NTIS, file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                                      "_", data, "HUGO_ntis.RData"))
          #bar chart
          p <- myplot(as.data.frame(HUGO.NTIS)) + geom_bar(stat="count")
          
        } else {
          
          #with separation of annotation types (NM, NR)
          #combine NM and NR tables differentiating between the two by adding column indicating if NM or NR
          HUGO.NTIS.NM <- cbind(Grp=rep("NM", nrow(HUGO.NTIS.list$NM)), HUGO.NTIS.list$NM)
          HUGO.NTIS.NR <- cbind(Grp=rep("NR", nrow(HUGO.NTIS.list$NR)), HUGO.NTIS.list$NR)
          HUGO.NTIS.NMNR <- rbind(HUGO.NTIS.NM, HUGO.NTIS.NR)
          
          save(HUGO.NTIS.NMNR, file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                                           "_", data, "HUGO_ntis.RData"))
          #bar chart
          p <- myplot(as.data.frame(HUGO.NTIS.NMNR)) + geom_bar(stat="count", aes(fill=Grp))
        }
        
        #for saving as pdf
        print(p)
        dev.off()
        
        ggsave(file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                           "_", data, "HUGO_ntis.png"))
        
      } #end foreach data for plotting
      
    } #end foreach chr
  } #end foreach gc
  
} else { #else statement for plotOnly
  
  #load objects
  foreach(gc=gc.v, .inorder=FALSE) %op% {
    foreach(chr=chr.v, .inorder=FALSE) %op% {
      foreach (data=c("ALL", "NMNR"), .inorder=FALSE) %op% {
        
        pdf(file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                        "_", data, "HUGO_ntis.pdf"), height=8, width=8)
        
        if(data=="ALL"){
        
          load(paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                      "_", data, "HUGO_ntis.Rdata"))
          p <- myplot(as.data.frame(HUGO.NTIS)) + geom_bar(stat="count")
          
        } else {
          
          load(paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                      "_", data, "HUGO_ntis.Rdata"))
          p <- myplot(as.data.frame(HUGO.NTIS.NMNR)) + geom_bar(stat="count", aes(fill=Grp))
          
        }
        
        #for saving as pdf
        print(p)
        dev.off()
        
        ggsave(file=paste0(output.dir, "/chr", chr, "_min", gc, "Mb_", genome.ver, 
                           "_", data, "HUGO_ntis.png"))
        
      } #end foreach data for plotting
    }
  }
  
} #end else statement of plotOnly

#rm(list=ls())
