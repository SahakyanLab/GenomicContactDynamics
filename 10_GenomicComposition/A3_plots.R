################################################################################
# Make volcano plots and sequence logos (all enriched/depleted kmers and also
# top5% enriched/depleted kmers based on alpha)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    # UNIQBIN.DF directory
    data.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/10_GenomicComposition"
  } else if(whorunsit == "LiezelCluster"){
    # UNIQBIN.DF directory
    data.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/9_GenomicComposition"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(wk.dir, "/out_compare")
out.dir = paste0(wk.dir, "/out_plots")
### OTHER SETTINGS #############################################################
gcb = "min05Mb"
kmer.len = "7"
# CHANGE x limits of volcano plot appropriately
cp.v = 1:21
kmerDistVal = "Mean"
SL.method = "prob" # "bit"
HM = FALSE
VPSL = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggpubr)
library(ggseqlogo) # geom_logo
library(ggrepel) # geom_text_repel
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/makeSeqLogo.R"))

transformFUN <- function(x) formatC(x, format="e", digits=1)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
#for(kmerDistVal in kmerDistVal.v){
  
  # Load KMERCP list
  load(file=paste0(data.dir, "/KMERCP", kmer.len, "_", gcb, "_", kmerDistVal,
                   ".RData")) 
  
  # NaNs will be introduced
  centrR.mx <- KMERCP$CENTR[,-1]/KMERCP$CENTR[,"HiCAll"] 
  kmer.v <- rownames(centrR.mx)
  # Infinite log10(pval) <- significant
  KMERCP$NEGLOG10PVAL[is.infinite(KMERCP$NEGLOG10PVAL)] <- ceiling( max(KMERCP$NEGLOG10PVAL[is.finite(KMERCP$NEGLOG10PVAL)], 
                                                                        na.rm=TRUE) ) + 1L
  
  # Check for not finite values including Nas, NaNs, ±Inf
  notFin <- unlist( lapply(X=KMERCP, FUN=function(x)sum(!is.finite(x)) )
                    )
  if( any(notFin)!=0L ){ stop("Checkpoint 1.") }
  rm(notFin); gc()
  
  cp.v <- as.character(cp.v)
  
  #---------------------------------------
  # Volcano plots and Sequence logos
  #---------------------------------------
  if(VPSL==TRUE){
    
    col.v <- c("gray50", "darkred", "darkblue", "darkgreen")
    names(col.v) <- c("n.s.", "sig.greater", "sig.less", "sig.equal")
    label.v <- c( "n.s.", expression( bold("greater") ),
                  expression( bold("less") ),
                  expression( bold("equal")) )
    names(label.v) <- c("n.s.", "sig.greater", "sig.less", "sig.equal")
    
    p.lst <- list()
    VOLCANO <- list()
    for(cp in cp.v){
      
      #---------------------------------------
      # Volcano plots
      #---------------------------------------
      VOLCANO[[cp]] <- data.frame( kmer=kmer.v,
                                   log2centrR=log2(centrR.mx[,cp]),
                                   neglog10pval=KMERCP$NEGLOG10PVAL[,cp],
                                   row.names=NULL, stringsAsFactors=FALSE
      )
      # NaN = 0/0, possible for both Median and Mean cases
      # Inf = x>0/0, will only happen for the Median case; impossible
      # for the mean case because sample can never be >0 if control=0
      # Remove kmers with log2centrR=NaN/±Inf 
      VOLCANO[[cp]] <- VOLCANO[[cp]][!is.nan(VOLCANO[[cp]][,"log2centrR"]),]
      #VOLCANO[[cp]] <- VOLCANO[[cp]][is.finite(VOLCANO[[cp]][,"log2centrR"]),]
      
      # -log10(0.05) ~ 1.301; The lower the alpha (meaning the more significant 
      # it is), the higher the neglog10pval
      # Positive log2centrR means samp.centr > cont.centr
      VOLCANO[[cp]] <- within(data=VOLCANO[[cp]], {
        group <- NA
        group[neglog10pval < 1.301] <- "n.s."
        group[neglog10pval >= 1.301 & log2centrR > 0L] <- "sig.greater"
        group[neglog10pval >= 1.301 & log2centrR < 0L] <- "sig.less"
        group[neglog10pval >= 1.301 & log2centrR==0L] <- "sig.equal"
      })
      
      # Check for NAs
      if( sum( is.na(VOLCANO[[cp]]) )!=0L ){ stop("Checkpoint 2.") }
      
      group.v <- levels( as.factor(VOLCANO[[cp]][,"group"]) )
      
      # Only label and make sequence logos, when there are kmers significantly
      # enriched and depleted
      log <- ifelse("sig.greater"%in%group.v & "sig.less"%in%group.v, TRUE, FALSE)
      if(log){
        VOLCANO[[cp]] <- VOLCANO[[cp]][order(VOLCANO[[cp]][,"group"], 
                                             VOLCANO[[cp]][,"neglog10pval"],
                                             decreasing=TRUE),]
        ind.less <- which(VOLCANO[[cp]][,"group"]=="sig.less")
        ind.less.len <- length(ind.less)
        ind.greater <- which(VOLCANO[[cp]][,"group"]=="sig.greater")
        ind.greater.len <- length(ind.greater)
        # Indices of kmers to label
        kmerlabind <- c(ind.less[1:5], ind.greater[1:5])
        
        VOLCANO[[cp]] <- within(data=VOLCANO[[cp]], {
          putlabel <- NA
          putlabel[kmer%in%VOLCANO[[cp]][kmerlabind,"kmer"] ] <- 1
        })
        
      }
      
      #max.y <- max(VOLCANO[[cp]][,"neglog10pval"])
      #max.x <- max(VOLCANO[[cp]][,"log2centrR"])
      
      # Number of significantly enriched/depleted kmers out of 16384 (4^7)
      totkmersSig <- sum(VOLCANO[[cp]][,"group"]!="n.s.")  
      id <- bquote(.(gcb)~"_"~.(kmer.len)~"mer_"~.(kmerDistVal)~"R_cp="~.(cp)~"VsHiCAll_"~alpha~"<0.05_"~.(totkmersSig)~"sigkmr")
      
      p.lst[[cp]] <- ggplot(data=VOLCANO[[cp]],
                            aes(x=log2centrR, y=neglog10pval)) +
        geom_point( size=1 ,aes(colour=as.factor(group)) ) +
        # CHANGE
        #scale_x_continuous(limits=c( floor( min(VOLCANO[[cp]][,"log2centrR"],
        #                                        na.rm=TRUE) ), 
        #                             ceiling(max.x)
        #) ) +
        scale_x_continuous(limits=c(-1.5,1.5)) +
        scale_y_continuous(labels=transformFUN) + 
        scale_colour_manual(values=col.v[ names(col.v)%in%group.v ], 
                            labels=label.v[ names(label.v)%in%group.v ]
        ) + 
        #guides( colour=guide_legend(override.aes=list(size=5)) ) +
        labs(title=id , 
             x=bquote(bold("log"["2"]~"("~.(kmerDistVal)~" FC)")), 
             y=bquote(bold("-log"["10"]~"(p-value)")),
             colour=""
        ) +
        #theme(legend.text=element_text(size=20, face="bold"),
        #      legend.title=element_text(size=25, face="bold"),
        #      legend.position="top",
        #      aspect.ratio=1) + 
        bgr2 +
        labs(title=NULL, x=NULL, y=NULL) + 
        theme(legend.position="none", aspect.ratio=1,
              axis.text.x=element_blank(),
              axis.text.y=element_text(size=15))
      
      if(log){
        
        p.lst[[cp]] <- p.lst[[cp]] + 
          geom_text_repel(data=subset(VOLCANO[[cp]], putlabel==1),
                          aes(label=kmer), size=3,
                          box.padding=unit(0.35, "lines"),
                          point.padding=unit(0.3, "lines"),  
                          max.overlaps=20 
          )
        
      }
      
      out.name <- paste0(out.dir, "/KMERCP", kmer.len, "_", gcb, "_cp", 
                         cp, "_", kmerDistVal)
      # Save individual volcano plots
      ggsave(filename=paste0(out.name, "_volcp.pdf"), 
             width=10, height=10, plot=p.lst[[cp]])
      
      #save(VOLCANO, file=paste0(out.name, "_volcp.RData"))
      
      #---------------------------------------
      # Sequence logos
      #---------------------------------------
      # Only make sequence logos when there are kmers significantly enriched and depleted
      if(log){
        
        DNAseq <- list()
        
        lab <- c(paste0( "Depleted_", round(100*ind.less.len/totkmersSig, digits=2), "%" ),
                 paste0( "Enriched_", round(100*ind.greater.len/totkmersSig, digits=2), "%" )
        )
        DNAseq[[lab[1]]] <- as.character(VOLCANO[[cp]][ind.less,"kmer"])
        DNAseq[[lab[2]]] <- as.character(VOLCANO[[cp]][ind.greater,"kmer"])
        
        p.seqL <- makeSeqLogo(list=DNAseq, ncol=2L, title=id, method=SL.method)
        
        ggsave(filename=paste0(out.name, "_seqLogo.pdf"), 
               width=10, height=8, plot=p.seqL)
        
        # Top 5% enriched/depleted
        DNAseq[[lab[1]]] <- DNAseq[[lab[1]]][1:trunc(ind.less.len*0.05)] 
        DNAseq[[lab[2]]] <- DNAseq[[lab[2]]][1:trunc(ind.greater.len*0.05)]
        p.seqL <- makeSeqLogo(list=DNAseq, ncol=2L, title=bquote(.(id)~"_top5"~alpha), method=SL.method)
        
        ggsave(filename=paste0(out.name, "_top5_seqLogo.pdf"), 
               width=10, height=8, plot=p.seqL)
        
        rm(DNAseq, p.seqL, ind.less, ind.less.len, ind.greater, 
           ind.greater.len, kmerlabind); gc()
        
      }
      
      rm(group.v, log, out.name, totkmersSig); gc()
      
      print(paste0("cp=", cp, " plot done!"), quote=FALSE)
      
    } # cp.v for loop end
    
    # Save volcano plots as a group
    p.lst <- ggarrange(plotlist=p.lst, nrow=3, ncol=7,
                       common.legend=FALSE
                       #, legend="top"
                       )
    #ggexport(p.lst, width=70, height=30,
    ggexport(p.lst, width=35, height=15,
             filename=paste0(out.dir, "/KMERCP", kmer.len, "_", gcb, 
                             "_", kmerDistVal, "_volcp.pdf"))
    save(VOLCANO, file=paste0(out.dir, "/KMERCP", kmer.len, "_", gcb, 
                              "_", kmerDistVal, "_volcp.RData"))
    
  }
  
  #---------------------------------------
  # Heatmaps
  #---------------------------------------
  if(HM==TRUE){
    
    coul <- c( colorRampPalette(c("darkblue", "blue", "red", "darkred"))(n=299),
               "darkseagreen3" )
    
    # Cluster kmers using central value ratios
    pdf(file=paste0(out.dir, "/KMERCP", kmer.len, "_", gcb, 
                    "_", kmerDistVal, "ratioHM.pdf"), width=10, height=8)
    HM <- heatmap.2(x=centrR.mx, Rowv=TRUE, Colv=FALSE, dendrogram="row",
                    scale="none", trace="none", na.rm=TRUE,
                    xlab=expression( bold("c"["p"]) ), 
                    ylab=expression( bold("7-mer") ), 
                    margins=c(4,5.5), col=coul, key=TRUE, key.title=NA, 
                    distfun = function(x) dist(x, method = "euclidean") )
    save(HM, file=paste0(out.dir, "/KMERCP", kmer.len, "_", gcb, 
                         "_", kmerDistVal, "ratioHM.RData"))
    dev.off()
    
    # Cluster kmers using p-values
    pdf(file=paste0(out.dir, "/KMERCP", kmer.len, "_", gcb, 
                    "pvalueHM.pdf"), width=10, height=8)
    HM <- heatmap.2(x=KMERCP$NEGLOG10PVAL, Rowv=TRUE, Colv=FALSE, dendrogram="row",
                    scale="none", trace="none", na.rm=TRUE,
                    xlab=expression( bold("c"["p"]) ), 
                    ylab=expression( bold("7-mer") ), 
                    margins=c(4,5.5), col=coul, key=TRUE, key.title=NA, 
                    distfun = function(x) dist(x, method = "euclidean") )
    save(HM, file=paste0(out.dir, "/KMERCP", kmer.len, "_", gcb, 
                         "pvalueHM.RData"))
    dev.off()
    
  }
  
#} # kmerDistVal.v for loop end 

# rm(list=ls()); gc()

