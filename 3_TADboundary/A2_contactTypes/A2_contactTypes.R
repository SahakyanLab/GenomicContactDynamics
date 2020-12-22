################################################################################
# Using the TAD boundary information, determine per cell/tissue the percentage 
# of the following:
# a. contacts between TAD boundaries
# b. countacts between TAD boundary and non-TAD boundary
# c. contacts between non-TAD boundaries (can be inter=TAD or intra-TAD)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/3_TADboundary"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/3_TADboundary"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# Bed files directory
TADb.dir = paste0(wk.dir, "/Schmitt2016_TADboundary")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_contactTypes")
### OTHER SETTINGS #############################################################
#ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
#              "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC"))
ct = "FC" #"REPLACE"
gcb.v = "min2Mb"
HiC.res = 4e4L
# Cp
nCPU = 5L #~26G
plotOnly = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
library(reshape)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Read bed file
bound <- fread(file=paste0(TADb.dir, "/ct_", ct, "_foi_TADboundary_desc_Schmitt"),
               select=c(1,3), col.names=c("chr", "bin"), 
               stringsAsFactors=FALSE, data.table=FALSE, header=FALSE)

# Convert TAD boundaries to bins
bound$bin <- bound$bin/HiC.res
chr.v <- unique(bound$chr)

for(gcb in gcb.v){
  
  p.lst <- list()
  
  for(ct in ct.v){
    
    if(plotOnly==FALSE){
      
      # Initialize output matrix
      CONTTYPE.MX <- matrix(data=0, nrow=21, ncol=5, 
                            dimnames=list(1:21, 
                                          c("b-b", "nb-b", "nb-nb(inter)", "nb-nb(intra)", "total"))
      )
      
      for(chr in chr.v){
        
        # Load PERSIST.MX
        load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
        # Only contacts in given cell/tissue
        log.ct <- PERSIST.MX$hits[[ct]]!=0
        ij.mx <- cbind(i=PERSIST.MX$hits[log.ct,"i"], 
                       j=PERSIST.MX$hits[log.ct,"j"],
                       Cp=PERSIST.MX$ntis[log.ct])
        rm(PERSIST.MX, log.ct); gc()
        
        cp.v <- sort(unique(ij.mx[,"Cp"], decreasing=FALSE))
        
        toExport <- c("ij.mx", "cp.v", "bound")
        
        #### PARALLEL EXECUTION #########
        temp <- foreach(itr=isplitVector(1:21, chunks=nCPU), 
                        .inorder=TRUE, .combine="rbind",
                        .export=toExport, .noexport=ls()[!ls()%in%toExport]
                        
        ) %op% {
          
          chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
            
            if(!i%in%cp.v){
              return( c(0,0,0,0,0) )
            } else {
              
              log.cp <- ij.mx[,"Cp"]==i
              log.bound <- bound$chr==chr
              log.ibound <- ij.mx[,"i"]%in%bound[log.bound, "bin"]
              log.jbound <- ij.mx[,"j"]%in%bound[log.bound, "bin"]
              
              ij.bound <- cbind(i=ij.mx[log.cp, "i"]%in%bound[log.bound, "bin"],
                                j=ij.mx[log.cp, "j"]%in%bound[log.bound, "bin"])
              
              # "TADb-TADb"
              a <- sum(log.cp&log.ibound&log.jbound)
              
              # "nTADb-TADb"
              b <- sum( apply(X=cbind(log.ibound[log.cp], log.jbound[log.cp]), 
                              MARGIN=1, FUN=function(x){
                                # Order not important
                                setequal(x, c(TRUE, FALSE))
                              }) )
              
              # "nTADb-nTADb"
              log.nbnb <- log.cp & !log.ibound & !log.jbound
              
              # inter "nTADb-nTADb"
              c1 <- sum(
                apply(X=ij.mx[log.nbnb, c("i","j")],
                      MARGIN=1, FUN=function(x){
                        return(any( 
                          ( (x[1]+1):(x[2]-1) )%in%bound[log.bound, "bin"] 
                        ))
                      })
              )
              
              # intra "nTADb-nTADb"
              c2 <- sum(
                
                apply(X=ij.mx[log.nbnb, c("i","j")],
                      MARGIN=1, FUN=function(x){
                        return(!any( 
                          ((x[1]+1):(x[2]-1))%in%bound[log.bound, "bin"] 
                        ))
                        
                      })
                
              )
              
              tot <- sum(log.cp)
              out <- c(a, b, c1, c2, tot)
              if( sum(out[-5])!=tot | tot!=nrow(ij.bound) ){
                stop("Checkpoint 2.")
              }
              
              #print(paste0("Cp=", i, " done!"), quote=FALSE)
              return(out)
              
            }
          }) # chunk sapply end
          
          return( do.call("rbind", chunk) )
          
        }
        ### END OF PARALLEL EXECUTION ###
        
        CONTTYPE.MX <- CONTTYPE.MX + temp
        
        rm(temp, ij.mx, cp.v); gc()
        
        print(chr, quote=FALSE)
        
      } # chr.v for loop end
      
      CONTTYPE.MX <- (CONTTYPE.MX[,-5]/CONTTYPE.MX[,5])*100
      
      save(CONTTYPE.MX, 
           #     file=paste0(out.dir, "/chrALL_", gcb, "_", ct, "_contactType.RData"))
           file=paste0(out.dir, "/", ct, "_", gcb, "_contactType.RData"))
      
      # Check if percentages add to 100
      if( any( round(rowSums(CONTTYPE.MX), digits=0)!=100L ) ){
        stop("Checkpoint 2.")
      }
      
    } else {
      
      # Load CONTTYPE.MX
      #load(file=paste0(out.dir, "/chrALL", gcb, "_", ct, "_contactType.RData"))
      load(file=paste0(out.dir, "/", ct, "_", gcb, "_contactType.RData"))
      
    }
    
    # Plot
    df <- melt(CONTTYPE.MX)
    colnames(df) <- c("Cp", "type", "value")
    id <- paste0("chrALL_", gcb, "_", ct)
    
    coul <- colorRampPalette( rev(brewer.pal(11, "Spectral")) )(21)
    p.lst[[ct]] <- ggplot(data=df, aes(x=type, y=value, group=Cp)) +
      geom_point(aes(colour=factor(Cp)), size=8) +
      scale_y_continuous(breaks=c(25,50,75)) +
      scale_x_discrete(breaks=c("b-b", "nb-b", "nb-nb(inter)", "nb-nb(intra)"),
                       labels=c("b-b", "nb-b", "nb-nb\ninter", "nb-nb\nintra")) +
      scale_colour_manual(values=coul) +
      guides( colour=guide_legend( title=expression(bold("C"["p"]) ),
                                   ncol=1) 
      ) + 
      labs( title=id, y=expression(bold("% Contacts")), x=NULL ) +
      bgr2 +
      theme(axis.text.x = element_text(face="bold", size=20, 
                                       angle=360, colour="black")
      ) #+
      #labs(title=NULL, y=NULL, x=NULL) +
      #theme(legend.position="none") + 
      #theme(axis.text.x=element_blank(),
      #      axis.text.y=element_blank())
    
    ggsave(filename=paste0(out.dir, "/", id, "_contactType.pdf"),
           units="in", height=10, width=10, plot=p.lst[[ct]])
    
    print(ct, quote=FALSE)
    
    rm(CONTTYPE.MX); gc()
    
  } #ct.v for loop end
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=3, ncol=7,
                     legend=NULL)
  ggexport(p.arr, height=30, width=70,
  #ggexport(p.arr, height=15, width=35,
           filename=paste0(out.dir, "/", gcb, "_contactType_Schmitt.pdf" ))
  
} # gcb.v for loop end

# rm(list=ls()); gc()






  
  

