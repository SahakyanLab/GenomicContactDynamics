################################################################################
# Determine radial position of contacting bins based on Chrom3D structure
# deva, R/3.6.0-newgcc, gcc/4.9.2
# athena, R/3.6.1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/6_Location_Chrom3D"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/project/sahakyanlab/ltamon/DPhil/lib"
    wk.dir = "/project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/6_Location_Chrom3D"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "IMR90_LMNB1_GSE49341_hg19" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
data.dir = out.dir = paste0(wk.dir, "/out_ContactRadDist/", model.id)
### OTHER SETTINGS #############################################################
ploidy = "haploid"
gcb.v = "min2Mb"
chr.v = "chrALL" #paste("chr", c(1:22, "X"), sep="")
# If doParCHR = FALSE, parallel execution with gcb.v
doParCHR = FALSE 
nCPU = 1L

boxPlot = FALSE
densPlot = TRUE
corPlot = FALSE
 ggscattr = FALSE
 ntis.corPlot = 1:21
 ntis.lab = "PScoreALL"
corCoefPlot = FALSE
plotOnlyCoeff = FALSE
cor.method = "spearman"
## Note
# Adjust to x limits of density plot to make models comparable
# Currently c(0,6) for both FC and ESC cell lines
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
out.name <- paste0(model.id, "_", ploidy) 
library(ggplot2)
if(corPlot==TRUE & ggscattr==TRUE){ library(ggpubr) }
library(RColorBrewer)
library(foreach)
library(itertools) #isplitVector
library(doParallel) 
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/GG_bgr.R"))

# Divided by 2 because each contact has two radial position values
binN <- function(x){ return(c(y=-0.3, label = length(x)/2)) }

myboxplot <- function( df=df,
                       uniqueNtis=unique(IJ.FEAT.MX[,"ntis"]),
                       title=plottitle
                      ){
  
  coul = colorRampPalette(rev(coul0))(length(uniqueNtis))
  
  ggplot(data=as.data.frame(df),
         aes(x=ntis, y=value, group=ntis)) +
    stat_boxplot(geom="errorbar", width = 0.3, coef = 1.5) + 
    geom_boxplot(fill="#FDC776", na.rm=FALSE) +
    stat_summary(fun.data=binN, geom="text", size=1) + 
    scale_x_discrete(limits=sort(as.numeric(uniqueNtis), decreasing=FALSE),
                     drop=FALSE) +
    #scale_fill_manual(values=coul) +
    guides(fill="none") + 
    labs(title=title, x=expression("c"["p"]), y=expression("r"["ij"])) +
    bgr2 
}

mydensplot <- function( df=df,
                        uniqueNtis=unique(IJ.FEAT.MX[,"ntis"]),
                        title=plottitle
                       ){
  
  coul = colorRampPalette(rev(coul0))(length(uniqueNtis))
  
  ggplot(data=as.data.frame(df),
         aes(value, colour=factor(ntis))) +
    geom_density(na.rm=FALSE) +
    # Adjust to make density plot comparable among models
    scale_x_continuous(limits=c(0,6), breaks=0:6) +
    #scale_y_continuous(limits=c(0,0.8), breaks=seq(0,0.8,0.2)) +
    scale_colour_manual(values=coul) +
    guides(colour=guide_legend(ncol=1)) + 
    labs(title=title, x="r", y="Density", colour=expression("c"["p"])) +
    bgr2 +
    theme( legend.text=element_text(size=40, face="bold") )
}

mycorplot <- function( df=df, title=plottitle ){
  
  RadDist.lims <- range(c(df[,"iRadDist"], df[,"jRadDist"]), na.rm=T)
  RadDist.lims[1] <- floor(RadDist.lims[1])
  RadDist.lims[2] <- ceiling(RadDist.lims[2])
  
  axis.brks <- RadDist.lims[1]:RadDist.lims[2]
  
  if(ggscattr==FALSE){
    ggplot(data=as.data.frame(df), aes(x=iRadDist, y=jRadDist)) +
      geom_point(color="gray59", shape=1) +
      stat_smooth(geom="smooth", method="lm", formula=y~x, se=TRUE, 
                  n=80, span=0.75, level=0.95, aes(colour="darkred"), size=3) +
      #scale_y_continuous(limits=c(-1, NA)) +
      scale_y_continuous(breaks=axis.brks) + 
      scale_x_continuous(breaks=axis.brks) + 
      guides(colour="none") + 
      labs( title=title, x=expression("r"["i"]), y=expression("r"["j"]) ) + 
      bgr3 + 
      facet_wrap( ~ ntis, ncol=3) + 
      theme(strip.text.x = element_text(size=40, angle=360, face="bold")
            #strip.background = element_rect(colour="gray22", fill="white")
            )
  } else {
    # ggscatter way
    ggscatter(as.data.frame(df), x="iRadDist", y="jRadDist", 
              color = "gray59", shape=1, 
              add="reg.line", cor.coef = FALSE, cor.method = "pearson",
              conf.int=TRUE, conf.int.level=0.95, point=TRUE, size=3,
              add.params=list(color = "darkred", fill = "red", size=3),
              xlab=expression("r"["i"]), ylab=expression("r"["j"]),
              title=title, facet.by="ntis") + 
      stat_cor(label.y=0, size=10, colour="darkred", method="pearson") +
      bgr2 
  }
  
}

mycorCoefPlot <- function(df=df, title=plottitle, cor.method){
  ggplot(data=df, aes(x=ntis, y=Coef, label=pvalue)) +
    #geom_line(size=2.5) +
    geom_point(size=6, colour="darkred") +
    scale_x_continuous(breaks=ntis.corPlot) +
    scale_y_continuous(limits=c(0,1)) +
    labs(title=title, x=expression("c"["p"]), 
         y=paste0(cor.method, " coefficient")) +
    bgr2 
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
coul0 <- brewer.pal(11, "Spectral")

out.name <- paste0(out.name, "_ContactRadDist")

toExport <- c("chr.v", "gcb.v", "out.dir", "out.name", "out.name", "ntis.corPlot", 
              "ntis.lab", "boxPlot", "densPlot", "corPlot", "ggscattr", 
              "corCoefPlot", "plotOnlyCoeff", "doParCHR")

if(doParCHR==TRUE){
  foreach.v.len <- length(chr.v)
  for.v <- gcb.v
} else {
  foreach.v.len <- length(gcb.v)
  for.v <- chr.v
}

for(v in for.v){
  
  toExport <- unique(c(toExport, "v"))
  
  foreach( itr=isplitVector(1:foreach.v.len, chunks=nCPU), 
           .inorder=FALSE, .export=toExport, 
           .noexport=ls()[!ls()%in%toExport]
           
  ) %op% {
    
    for(i in itr){
      
      if(doParCHR==TRUE){
        chr <- chr.v[i]
        gcb <- v
      } else {
        gcb <- gcb.v[i]
        chr <- v
      }
      
      # Load IJ.FEAT.MX
      if( any(c(boxPlot, densPlot, corPlot, corCoefPlot)) ){
        # Load IJ.FEAT.MX
        load(file=paste0(data.dir, "/", chr, "_", gcb, "_", out.name, 
                         ".RData"))
        ntis.uniq <- unique(IJ.FEAT.MX[,"ntis"])
        ntis.uniq.len <- length(ntis.uniq)
      }
    
      plottitle <- paste0(chr, "_", gcb, "_", out.name)
      
      # Make boxplot radDistVsCp
      if(boxPlot==TRUE){
        # Contacts are unique so rbind is fine
        mx <- rbind(IJ.FEAT.MX[,c("ntis", "iRadDist")],
                    IJ.FEAT.MX[,c("ntis", "jRadDist")])
        dimnames(mx) <- list(NULL, c("ntis", "value"))
        # na.omit() to only include contacts based on setAgeCountPerBin
        myboxplot( df=na.omit(mx), uniqueNtis=ntis.uniq, title=plottitle )
        ggsave(filename=paste0(out.dir, "/", chr, "_", gcb, "_", 
                               out.name, "_bp.jpeg"), 
               units="in", width=13, height=13)
        rm(mx)
      }
      
      if(densPlot==TRUE){
        mx <- rbind(IJ.FEAT.MX[,c("ntis", "iRadDist")],
                    IJ.FEAT.MX[,c("ntis", "jRadDist")])
        dimnames(mx) <- list(NULL, c("ntis", "value"))
        # na.omit to only include contacts based on setAgeCountPerBin
        mydensplot( df=na.omit(mx), uniqueNtis=ntis.uniq, title=plottitle )
        ggsave(filename=paste0(out.dir, "/", chr, "_", gcb, "_", 
                               out.name, "_dp.jpeg"), 
               units="in", width=10, height=10)
        rm(mx)
      }
      
      # Make correlation plot radial Positions of c(i,j)
      if(corPlot==TRUE){
        
        test1 <- IJ.FEAT.MX[,"ntis"]%in%ntis.corPlot
        # na.omit to only include contacts based on setAgeCountPerBin
        cp <- mycorplot( df=na.omit(IJ.FEAT.MX[test1,c("ntis", "iRadDist", "jRadDist")]),
                         title=plottitle) 
        ggsave(file=paste0(out.dir, "/", chr, "_", gcb, "_", ntis.lab, "_", 
                           out.name, "_cp.jpeg"), units="in", width=20, height=35, plot=cp)
        
        test2 <- IJ.FEAT.MX[,"ntis"]%in%c(1, 5, 9, 13, 17, 21)
        cp1 <- mycorplot( df=na.omit(IJ.FEAT.MX[test2,c("ntis", "iRadDist", "jRadDist")]),
                          title=plottitle) 
        ggsave(file=paste0(out.dir, "/", chr, "_", gcb, "_2by3_", out.name, "_cp.jpeg"), 
               units="in", width=20, height=15, plot=cp1)
        
      }
      
      if(corCoefPlot==TRUE){
        
        if(plotOnlyCoeff==FALSE){
          
          corr <- sapply(X=ntis.corPlot, simplify=FALSE, FUN=function(ntis){
            test1 <- IJ.FEAT.MX[,"ntis"]==ntis
            corr <- cor.test(IJ.FEAT.MX[test1,"iRadDist"], 
                             IJ.FEAT.MX[test1, "jRadDist"],
                             method=cor.method, exact=FALSE, conf.level=0.95)
            # Value lower than .Machine$double.eps (2.220446e-16) is returned as 0
            ifelse( corr$p.value==0, pval <- "<2.2e-16", 
                    pval <- formatC(x=corr$p.value, format="e", digits=1) )
            data.frame(ntis=ntis, Coef=corr$estimate, pvalue=pval, pvalue.actual=corr$p.value,
                       row.names=NULL, stringsAsFactors=FALSE)
          })
          
          df <- do.call(what="rbind", corr)
          save(df, file=paste0(out.dir, "/", chr, "_", gcb, "_", 
                               out.name, "_", cor.method, "_PCoeff.RData"))
        } else {
          load(file=paste0(out.dir, "/", chr, "_", gcb, "_", 
                           out.name, "_", cor.method, "_PCoeff.RData"))
        }
        
        cfp <- mycorCoefPlot(df=df, title=plottitle, cor.method=cor.method)
        ggsave(filename=paste0(out.dir, "/", chr, "_", gcb, "_", 
                               out.name, "_", cor.method, "_PCoeff.jpeg"), 
               units="in", width=13, height=13)
      }
      
      print(paste0(gcb, chr, ":", "ContactRadDist done!"), quote=FALSE)
      
    } # itr for loop end
    
  } 
} 

# rm(list=ls()); gc()




