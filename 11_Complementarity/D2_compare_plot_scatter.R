################################################################################
# Plot Cf vs. C|| (kmer and lign)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/11_Complementarity"
  } else if(whorunsit == "LiezelCluster"){
    #lib = "/t1-data/user/ltamon/DPhil/lib"
    #wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/pending/11_Constraints"
    lib = "/stopgap/sahakyanlab/ltamon/DPhil/lib"
    wk.dir = "/stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/pending/11_Constraints"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
Cs.form = "HiCNormCs" #"rawCs" | "HiCNormCs"
data.dir = paste0(wk.dir, "/out_compare_", Cs.form)
out.dir = paste0(wk.dir, "/out_compare_plot_scatter")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chrALL"
ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
              "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "LC", "FC"))
#ct.v = "Lu"
nCPU = 1
type = "Gfree" # "kmer" | "align" | "Gfree"
res = 300
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
x.col <- paste0("CII", type)
toExport <- c(Cs.form, data.dir, out.dir, gcb, chr, ct.v, type, x.col)

#### PARALLEL EXECUTION #########

ct.v.len <- length(ct.v)
P.LST <- foreach(itr=isplitVector(1:ct.v.len, chunks=nCPU), .inorder=TRUE, .combine="c",
                 .export=toExport, .noexport=ls()[!ls()%in%toExport]
      
) %op% {
  
  chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
    
    ct <- ct.v[i]
    id <- paste0(chr, "_", gcb, "_", ct)
    
    # Load CSCPCII.MX (actually a dataframe)
    load(file=paste0(data.dir, "/", id, "_CsCpCII.RData"))
    CSCPCII.MX$Cp <- as.numeric(as.character(CSCPCII.MX$Cp))
    
    len <- length(CSCPCII.MX[,1])
    
    # Include only contacts that have values for both metrics
    nonNA.TF <- !is.na(CSCPCII.MX[,x.col]) & !is.na(CSCPCII.MX$Cs)
    CSCPCII.MX <- CSCPCII.MX[nonNA.TF, c(x.col, "Cs", "Cp")]
    
    LR.TF <- !is.na(CSCPCII.MX$Cp)
    SR.TF <- is.na(CSCPCII.MX$Cp)
    nonNA.len <- length(LR.TF)
    CSCPCII.MX$Cp <- NULL
    
    num.ij <- c(tot.ij=len, nonNA.ij=sum(nonNA.TF), nonNALR.ij=sum(LR.TF), nonNASR.ij=sum(SR.TF))/len*100
    num.ij <- format(x=num.ij, scientific=FALSE, digits=4)
    num.ij <- paste0("percRelToTot.ij_", paste(paste(names(num.ij), num.ij, sep="_"), collapse="_"))
    rm(nonNA.TF)
    
    ijset <- rep(x="SR", times=nonNA.len)
    ijset[LR.TF] <- "LR"
    ijset <- factor(x=ijset, levels=c("LR", "SR"))
    
    p.id <- paste0(id, "_", Cs.form, "_", type)
    
    p <- ggplot(data=data.frame(CSCPCII.MX, ijset=ijset, stringsAsFactors=FALSE), 
                aes_string(x=x.col, y="Cs")) + 
      geom_point(colour="gray80", alpha=0.5, shape=1) +
      #geom_hex(bins=50) +
      geom_smooth(data=data.frame(CSCPCII.MX, ijset=factor("LR", levels=c("LR", "SR")), stringsAsFactors=FALSE), 
                  colour="#F39B7FFF", formula=y~x, method="lm", size=3, se=TRUE, n=80, level=0.95) +
      geom_smooth(data=data.frame(CSCPCII.MX, ijset=factor("SR", levels=c("LR", "SR")), stringsAsFactors=FALSE), 
                  colour="#F39B7FFF", formula=y~x, method="lm", size=3, se=TRUE, n=80, level=0.95) +
      geom_smooth(aes(colour=ijset), formula=y~x, method="lm", size=3, se=TRUE, n=80, level=0.95) +
      scale_colour_manual(values=c("#00A087FF", "#8491B4FF")) +
      guides(colour="none") + 
      labs(title=paste0(p.id, "_", num.ij), x=x.col, y="Cs") +
      bgr2 + 
      facet_grid(.~ijset)
      
    ggsave(filename=paste0(out.dir, "/", p.id, ".png"), width=30, height=15, 
           units="in", dpi=100, plot=p)
    
    # Calculate correlations
    ijset <- c("All", "SR", "LR")
    COR <- foreach(ijs=ijset, .combine="rbind"
                            
    ) %do% {
      
      if(ijs=="All"){
        ij.TF <- rep(x=TRUE, times=nonNA.len)
      } else {
        eval(parse(text=paste0( 'ij.TF <- ', ijs, '.TF' )))
      }
      
      eval(parse(text=paste0("lm.obj <- lm(Cs~", x.col, ", data=CSCPCII.MX[ij.TF,])")))
      coef.v <- setNames(object=coef(lm.obj)[c("(Intercept)", x.col)],
                         nm=c("intercept", "slope"))

      return(
        
        c(coef.v, Rsquared=summary(lm.obj)$r.squared,
          pearson=cor(x=CSCPCII.MX[ij.TF,x.col], y=CSCPCII.MX[ij.TF,"Cs"], method="pearson"),
          kendall=cor(x=CSCPCII.MX[ij.TF,x.col], y=CSCPCII.MX[ij.TF,"Cs"], method="kendall"),
          spearman=cor(x=CSCPCII.MX[ij.TF,x.col], y=CSCPCII.MX[ij.TF,"Cs"], method="spearman")) 
        
      )
      
    }
    
    rownames(COR) <- ijset
    
    rm(CSCPCII.MX, lm.obj, ijset)
    
    save(COR, file=paste0(out.dir, "/", p.id, "_cor.RData"))
    
    return(p)
    
  }) # ct.v sapply end
  
  return(chunk)
  
}

### END OF PARALLEL EXECUTION ###

P.LST <- lapply(X=P.LST, FUN=function(p){
  
  p <- p + labs(title=NULL, x=NULL, y=NULL)
  return(p)
  
})
p.arr <- ggarrange(plotlist=P.LST, nrow=3, ncol=7, legend=NULL)
ggexport(p.arr, height=5*3*res, width=10*7*res, res=res,
         filename=paste0(out.dir, "/", gcb, "_", chr, "_", Cs.form, "_", type, 
                         "_res", res, "_CsVsCII.png"))

# rm(list=ls()); gc()

