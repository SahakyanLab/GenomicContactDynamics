################################################################################
# Density plot of the nearest distance (log10 transformed) to chromosome ends 
# of the midPoints of unique bins per Cp (combined for all chr). 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/2_HiC_Human21_Expl_ext"
    persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_HiC_Human21_Expl_ext"
    persist.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = paste0(wk.dir, "/out_location")
# File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb.v = c("min2Mb", "min05Mb")
chr.v = paste("chr", c(22:1, "X"), sep="") 
# Bins per Cp per Chr
nCPU = 10L #~60G
HiC.res = 4e4L
plotOnly = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
library(ggplot2)
library(RColorBrewer)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(gcb in gcb.v){
  
  if(plotOnly==FALSE){
    
    # Chromosome length file
    chrLen.df <- fread(file=chrLenfile, colClasses=list(character=1, integer=2), 
                       stringsAsFactors=FALSE, data.table=FALSE)
    
    chr.v.len <- length(chr.v)
    
    DISTNEAREND <- replicate(n=21, simplify="list",
                             expr=matrix(, nrow=0, ncol=2))
    names(DISTNEAREND) <- 1:21
    
    getMidP <- (HiC.res-1)/2 
    
    for(c in 1:chr.v.len){
      
      chr <- chr.v[c]
      print(paste0(chr, "..."), quote=FALSE)
      
      maxPosChr <- as.integer(chrLen.df[chrLen.df[,1]==chr, 2])
      
      # Load PERSIST.MX
      load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      
      binsPerCp <- by(data=PERSIST.MX$hits[,c("i","j")],
                      INDICES=PERSIST.MX$ntis,
                      FUN=function(x){ unique( c(unique(x$i), unique(x$j)) ) }
      )
      rm(PERSIST.MX); gc()
      
      if(c==1){
        cp.v <- names(binsPerCp)
      }
      
      for(cp in cp.v){
        
        bins.v <- binsPerCp[[cp]]
        bins.v.len <- length(bins.v)
        
        #### PARALLEL EXECUTION #########
        temp <- foreach( itr=isplitVector(1:bins.v.len, chunks=nCPU),
                         .inorder=TRUE, .combine="c",
                         .export=c("bins.v", "HiC.res", "getMidP", "maxPosChr"), 
                         .noexport=ls()[!ls()%in%c("bins.v", "HiC.res", "getMidP", "maxPosChr")]
                         
        ) %op% {
          
          chunk <- sapply(X=itr, simplify=TRUE, FUN=function(i){
            # Exact midpoint of bins
            midPos <- bins.v[i]*HiC.res-getMidP
            distNE <- c(1-midPos, maxPosChr-midPos)
            distNE <- distNE[ which( abs(distNE)==min(abs(distNE)) ) ]
          })
          return(chunk)
        }
        ### END OF PARALLEL EXECUTION ###
        
        DISTNEAREND[[cp]] <- rbind( DISTNEAREND[[cp]],
                                    cbind(rep(as.numeric(cp)), 
                                          log10(abs(temp))) )
        colnames(DISTNEAREND[[cp]]) <-  c("cp", "log10dNE")
        rm(temp, bins.v); gc()
        
        print(paste0("cp=", cp), quote=FALSE)
        
      } # cp.v for loop end
      
      rm(binsPerCp); gc()
      print(paste0(chr, " done!"), quote=FALSE)
      
    } # chr.v.len for loop end
    
    DISTNEAREND <- do.call("rbind", DISTNEAREND)
    save(DISTNEAREND, file=paste0(out.dir, "/chrALL_distNE_", gcb, ".RData"))
    
  } else {
    
    load(file=paste0(out.dir, "/chrALL_distNE_", gcb, ".RData"))
    
  }
  
  coul <- colorRampPalette( rev( brewer.pal(11, "Spectral") ) )( length( unique(DISTNEAREND[,"cp"]) ))
  
  # Density plot
  p <- ggplot(data=as.data.frame(DISTNEAREND), aes(x=log10dNE)) +
    geom_density( position="identity", aes(colour=factor(cp)) ) + 
    scale_colour_manual(values=coul) + 
    labs( title=paste0("chrALL_distNE_", gcb), 
          x=bquote(bold("log"["10"]~"("~"d"^"NE"~")")),
          y=expression(bold(Density)), 
          colour=expression( bold( "c"["p"]) )
    ) +
    guides(colour=guide_legend(ncol=1)) +
    bgr2 
  
  ggsave(filename=paste0(out.dir, "/chrALL_distNE_", gcb, ".pdf"),
         units="in", width=10, height=10, plot=p)
  
}

# rm(list=ls())