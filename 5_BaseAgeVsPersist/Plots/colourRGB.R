#BaseAgeColRGBVsPersist - plotting the base age ColRGB
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
output.dir <- paste0(objective.dir, "/out_colourRGB") 

################################################################################
#Packages
library(ggplot2)
library(foreach)
library(RColorBrewer)

################################################################################
#Set values

#2(2Mb) or "05"(0-5 Mb), minimum gap accepted between contacting bins 
#gc.v <- c(2, "05")
gc <- 2
#chromosome
#chr.v <- c(1:22, "X")
chr <- "22"
#Column names in the BaseAgeVsPersist files
ColRGB <- "colourRGB"
Ntis <- "ntis"
plotOnly <- FALSE
#if TRUE (select colors)
colToPlot <- c("Reds", "Blues", "Greys") #or "Grey"

################################################################################   
#Functions

source(paste0(lib, "/ggLayers_persist.R"))

myplot <- function(dta) {
  coul = brewer.pal(9, i)
  coul = colorRampPalette(coul)( length( unique(dta[,"sumRGB"]) ) )
  
  ggplot(data=as.data.frame(dta), aes(x=dta[,Ntis], fill=dta[,"sumRGB"]))+
    geom_bar(stat="count")+
    scale_fill_manual(values=coul)+
    guides(fill=guide_legend(title="Shades"))+
    ggtitle(paste0("chr", chr, " _min", gc, "Mb"))+
    ylab("COUNT (Base)")+
    scale_x_discrete(name="COUNT (non-0 Tissues/Cell lines)", limits=c(1:21))+
    theme.persist  
}

################################################################################  
#foreach(gc=gc.v, .inorder=TRUE) %do% {
#  foreach(chr=chr.v, .inorder=TRUE) %do% {
    
    #names(colForPlot) <- c("Red", "Blue", "Grey")
    
    if(plotOnly==FALSE){
      load(paste0(baseAge.persist.dir, "/chr", chr, "_min", gc, 
                  "Mb_BaseAgeVsPersist.RData"))
      df <- BASEAGE.PERSIST$POOL.AND.MX
      
      if(is.null(df)){
        stop(paste0("No data for Chromosome ", chr, "."))
      }
      ColRGBPerNtis.lst <-  by( df[,ColRGB], df[,Ntis], 
                              function(x){
                                as.vector(x)
                              } )
      ColRGBPerNtis.lst <- lapply( as.list(names(ColRGBPerNtis.lst)), 
                             function(x){
                               perNtis <- ColRGBPerNtis.lst[[x]]
                               cllpseString <- paste(perNtis, collapse=";")
                               spltToRGBval <- strsplit(cllpseString, 
                                                        split=";")[[1]]
                               cbind(rep(x, length(spltToRGBval)), 
                                     spltToRGBval)
                                } )
      ColRGBPerNtis.df <- as.data.frame(do.call(rbind, ColRGBPerNtis.lst))
      #colnames(ColRGBPerNtis.df) <- c(Ntis, ColRGB)
      
      spltToR.G.B  <- sapply(ColRGBPerNtis.df[,2], 
                             function(y) strsplit(as.character(y), split=","))
      spltToR.G.B  <- do.call(rbind, spltToR.G.B)
      #colnames(spltToR.G.B) <- c("Red", "Green", "Blue")
      
      spltToR.G.B.sum <- apply(spltToR.G.B, MARGIN=1, 
                               function(x) sum (as.integer(x)))
      
      ind.addZeros <- which(nchar(spltToR.G.B.sum)!=3)
      if(length(ind.addZeros)!=0){
        addZeros <- spltToR.G.B.sum[ind.addZeros]
        len <- length(addZeros)
        
        withZeros <- c()
        for (i in 1:len){
          x <- addZeros[i]
          while(nchar(x)!=3){
            x <- paste0("0",x)
          }
          withZeros[i] <- x
        }
        
        spltToR.G.B.sum <- replace(spltToR.G.B.sum, ind.addZeros, withZeros)
      }
     
      RGBclassify <- apply(spltToR.G.B, MARGIN=1,
                           function(x){
                             if(length(unique(x))==1){
                               return("Greys")
                             } else {
                               if(max(as.numeric(x))==x[3]){
                                 return("Blues")
                               } else {
                                 ifelse(max(as.numeric(x))==x[1],  
                                        return("Reds"), 
                                        stop("Invalid shades present for 
                                             Chromosome ", chr, "."))
                               }
                             }
                           })
      
      BASEAGEcolRGB.PERSIST <- list()
      BASEAGEcolRGB.PERSIST <- lapply( list(Reds="Reds", Blues="Blues", 
                                            Greys="Greys"), 
        function(x){
          df <- cbind(ColRGBPerNtis.df[which(RGBclassify==x),],
                      sumRGB=spltToR.G.B.sum[which(RGBclassify==x)])
          rownames(df) <- NULL
          #df <- df[order(df[,Ntis], decreasing=FALSE),]
          #df <- df[order(df[,"sumRGB"], decreasing=TRUE),]
          df <- as.data.frame(df)
          df$sumRGB <- factor(df$sumRGB, levels=rev(levels(df$sumRGB)))
          df
        })
      
      save(BASEAGEcolRGB.PERSIST, 
           file=paste0(output.dir, "/chr", chr, "_min", gc, 
                       "Mb_BaseAgeColRGBVsPersist.RData"))
      
      for( i in names(BASEAGEcolRGB.PERSIST) ){
        pdf(file=paste0(output.dir, "/chr", chr, "_min", gc, 
                        "Mb_", i, "_BaseAgeColRGBVsPersist.pdf"),
            height=8, width=10)
        p <- myplot(BASEAGEcolRGB.PERSIST[[i]])
        print(p)
        dev.off()
      }
      
    } else {
      load(file=paste0(output.dir, "/chr", chr, "_min", gc, 
                       "Mb_BaseAgeColRGBVsPersist.RData"))
      
      for(i in colToPlot){
        pdf(file=paste0(output.dir, "/chr", chr, "_min", gc, 
                        "Mb_", i, "_BaseAgeColRGBVsPersist.pdf"),
            height=8, width=10)
        p <- myplot(BASEAGEcolRGB.PERSIST[[i]])
        print(p)
        dev.off()
      }
    
    }
#  }
#}

#end_time <- Sys.time()
#end_time-start_time   

#rm(list=ls())




