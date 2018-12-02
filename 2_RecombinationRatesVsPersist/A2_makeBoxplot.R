#Generate plots for recombination rates~persistent contacts 

################################################################################ 
#Set directories

#functions to source
lib <- "/Users/ltamon/ProjectBoard/lib"
#main directory for the task
objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/2_RecombinationRatesVsPersist"
#RecombinationRatesVsPersist files
recom.persist.dir <- paste0(objective.dir, "/out_A1_scr")
#dir.create(paste0(objective.dir, "/out_makeBoxplot"))
output.dir <- paste0(objective.dir, "/out_makeBoxplot")  

################################################################################  
#Packages
##CRAN
#install.packages("ggplot2")
library(ggplot2)

#install.packages("foreach") 
library(foreach)

################################################################################  
#Set values 

NCPU <- 1

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify 
#a contact, two points should be far enough to filter for contacts within a TAD
gc.v <- c(2, "05")

#write as a vector of chromosome number or X 
#no recombination data for Y and MT from Myers2005/Munsch2014
chr.v <- c("ALL") #c(1:22, "X", "ALL")

#rates in each bin were derived by taking the mean as well as median of the 
#individual rates in that bin
rb.v <- c("MEANrb", "MEDIANrb")

#column name of ntis (number of tissues) of ij pairs in RecombinationRatesVsPersist files
ntis <- "NTISij"

plotOnly <- FALSE #TRUE or FALSE
#if plotOnly <- TRUE, choose the plots you want to generate
#select.st should contain these chosen plots, st stands for the statistical method applied
#ij  - means statistical method using the mean or median of i and j bins contacting
#ijP - means statistical method using the the pooled values of i and j bins contacting
# select from these choices: c("MEANij", "MEDIANij", "SDEVij", "MINij", "MAXij", "DIFFij",
#                              "MEANijP", "MEDIANijP", "SDEVijP", "MINijP", "MAXijP")

select.st <- c("MEANij", "MEDIANijP")     
#change gc, chr, and rb accordingly

################################################################################   
#Functions

source(paste0(lib, "/UTL_NCPU.R"))
source(paste0(lib, "/ggLayers_persist.R"))

#give data points per box in boxplot
give.n <- function(x){ return(c(y=-0.5, label = length(x))) }

################################################################################ 
################################################################################   
#nested foreach to get data for each chromosome for each gap 
foreach(gc=gc.v, .inorder=FALSE) %op% {
  
  foreach(chr=chr.v, .inorder=FALSE) %op% {
    
    #load data
    load(paste0(recom.persist.dir, "/chr", chr, "_min", gc, "Mb_RecomVsPersist.RData"))
    
    foreach(rb=rb.v, .inorder=FALSE) %op% {
      
      #ggplot only takes a data frame
      df <- as.data.frame(RECOM.PERSIST$RECOM.PERSIST.MX[[rb]])
      #make sure data is not a list
      df <- as.data.frame(lapply(df, unlist))
      
      if(plotOnly==TRUE){
        
        #st.v should be set
        st.v <- select.st
    
      } else {
        
        if(rb=="MEANrb"){
          st.v <- c("MEANij", "MEDIANij", "SDEVij", "MINij", "MAXij", "DIFFij",
                    "MEANijP", "MEDIANijP", "SDEVijP", "MINijP", "MAXijP")
          
        } else {
          st.v <- c("MEANij", "MEDIANij", "SDEVij", "MINij", "MAXij", "DIFFij")
        }
      }
      
      #make RecombinationRatesVsPersist plots (per chromosome)
      foreach(st=st.v, .inorder=FALSE) %op% {
        
        #png(filename=paste0(recom.persist.dir, "filename"),
        #   width=800, height=640, res=64) #4x3 aspect ratio, 
        #default res: 72 pixels/inch
        #change res such that height/res=8-10 inches
        #when using png(), print ggplot using print()
        
        #subset data to be plotted
        data <- df[,c(ntis, st)]
        
        #remove rows with NAs
        NA.rows      <- apply(data, 1, function(x){any(is.na(x))})
        NA.rows.ind  <- which(NA.rows==TRUE)
        drop.rows.NA <- sum(NA.rows) #number of dropped rows (with NA)
        if(drop.rows.NA!=0){
          data <- data[-c(NA.rows.ind),]
        } 
        
        p <- ggplot(data=data, aes(x=data[,ntis], y=data[,st], group=data[,ntis])) +
          geom_boxplot(aes(fill=factor(ntis))) +
          guides(fill=FALSE) + 
          #use size 4 when using png (instead of ggsave)
          stat_summary(fun.data=give.n, geom="text", size=2) + 
          ggtitle(paste0("chr", chr, " _min", gc, "Mb", "_", st)) +
          #NEED TO FIX COLORS
          scale_x_discrete(name="COUNT (non-0/total Tissues)", 
                           limits=c(1:21), drop=FALSE) +
          #scale_x_continuous(name="COUNT (non-0/total Tissues)", 
          #                   breaks=unique(df[,"NTISij"]), 
          #                   labels=unique(df[,"NTISij"])) +
          scale_y_continuous(name=paste0(st, " RecombinationRate(cM/Mb)")) +
          theme.persist #from sourced "ggLayers_persist.R"
          
        if(drop.rows.NA!=0){
          p <- p + annotate("text", x=21, y=max(data[,st])+1.5, 
                            label=paste0("Dropped:", drop.rows.NA), size=3)
        } 
        
        if( st%in%c("MEANijP", "MEDIANijP", "SDEVijP", "MINijP", "MAXijP") ){
          
          pdf(file=paste0(output.dir, "/bp_min", gc, "Mb_", st, "_chr", chr, "_RecomVsPersist.pdf"),
              height=8, width=8)
          
          print(p)
          dev.off()
          
          ggsave(filename=paste0(output.dir, "/bp_min", gc, "Mb_", st, "_chr", chr, "_RecomVsPersist.png"),
                 plot=p)
        } else {
          
          pdf(file=paste0(output.dir, "/bp_min", gc, "Mb_", rb, "_", st,"_chr", chr, "_RecomVsPersist.png"),
              height=8, width=8)
          
          print(p)
          dev.off()
          
          ggsave(filename=paste0(output.dir, "/bp_min", gc, "Mb_", rb, "_", st, "_chr", chr, "_RecomVsPersist.png"),
                 plot=p)
        }

      } #end bracket for stat (st) foreach
    } #end bracket for type of rate per bin (rb) foreach
    
  } #end bracket for chromosome foreach
  
} #end bracket for gap foreach  

#rm(list=ls())




