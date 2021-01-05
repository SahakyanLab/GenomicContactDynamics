################################################################################
# Plot Cs vs. Cp, Cs vs. C||
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
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(wk.dir, "/out_compare_HiCNormCs")
out.dir = paste0(wk.dir, "/out_compare_plot/nolabel")
### OTHER SETTINGS #############################################################
ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
         "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "LC", "FC"))
gcb = "min2Mb"
chr = "chrALL"
ijset = "SR" # Contact set, All | LR (long-range) | SR (short-range)
bins=30
# Cs vs. Cp
CsCpPlot = FALSE
cuts=8
CsCIIPlot = TRUE
type.v = c("kmer", "align")
Cs.form = "HiCNorm" #c("raw", "log10", "HiCNorm")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(Hmisc)
library(hexbin)
library(viridis)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/makeHexbinggplot.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(CsCIIPlot){
  # Cs (raw and log10-transformed) vs. CII (kmer and align)
  y.v <- c(log10="log10(Cs)", raw="Cs", HiCNorm="Cs")
  y.lab <- list(raw=bquote(bold("raw "~"C"["s"])),
                HiCNorm=bquote(bold("HiCNorm "~"C"["s"])),
                log10=bquote(bold("log"["10"]~"C"["s"]))
                )
}
y.v <- y.v[Cs.form]

p.lst <- list()
p.kmer <- p.align <- list()

for(ct in ct.v){
  
  id <- paste0(chr, "_", gcb, "_", ct)
  load(file=paste0(data.dir, "/", id, "_CsCpCII.RData"))
  
  # Contact set, All | long-range | short-range
  if(ijset!="All"){
    
    if(ijset=="LR"){
      
      ij.TF <- !is.na(CSCPCII.MX[,"Cp"])
      print("Long-range contacts only...", quote=FALSE)
      
    } else if(ijset=="SR"){
      
      ij.TF <- is.na(CSCPCII.MX[,"Cp"])
      print("Short-range contacts only...", quote=FALSE)
      
    } else {
      stop("Wrong ijset argument!")
    }
    
  } else {
    
    ij.TF <- rep(TRUE, times=nrow(CSCPCII.MX))
    print("All contacts...", quote=FALSE)
    
  } 
  
  CSCPCII.MX <- CSCPCII.MX[ij.TF,]
  len <- nrow(CSCPCII.MX)
  tot.ij <- format(as.numeric(len), scientific=TRUE, digits=4)
  #---------------------------------------
  # Cs vs. Cp
  if(CsCpPlot){
    p <- makeHexbinggplot(xvar=CSCPCII.MX[,"Cp"], 
                          yvar=CSCPCII.MX[,"Cs"], 
                          bins=bins, 
                          cuts=cuts,
                          xlab=expression( bold("C"["p"]) ),
                          ylab=expression( bold("C"["s"]) ),
                          title=paste0(id, "_cuts", cuts, "_bins", bins, "_", 
                                       tot.ij, "ij"),
                          col=viridis(cuts)
    ) 
    ggsave(filename=paste0(out.dir, "/", id, "_", ijset, "_cuts", cuts, "_bins", 
                           bins, "_", Cs.form, "CsVsCp.pdf"),
           units="in", width=10, height=10, plot=p$hexplot)
    p.lst[[paste0(ct, "CsCp")]] <- p$hexplot; rm(p)
  }
  #---------------------------------------
  if(CsCIIPlot){
    
    for(type in type.v){
      nonNA.TF <- !is.na(CSCPCII.MX[[paste0("CII", type)]])
      nonNA.ij <- format(sum(nonNA.TF)/len*100, scientific=TRUE, 
                         digits=4)
      for(y.nme in names(y.v)){
        p.id <- paste0(ct, "_", type, "_cuts", cuts, "_bins", bins, "_", y.nme)
        eval(parse(text=paste0(
          "p <- ggplot(data=CSCPCII.MX[nonNA.TF,], aes(x=CII", type, ", y=", y.v[y.nme], "))"
        )))
        p <- p +
          geom_hex(bins=bins) +
          scale_fill_continuous(type="viridis") +
          #labs( y=y.lab[[y.nme]],
          #      x=bquote(bold( "C"["||"]~.(type) )),
          #      title=paste0(chr, "_", gcb, "_", p.id, "_", tot.ij, ijset, 
          #                   "ij_", nonNA.ij, "nonNACII") ) +
          bgr2 +
          #theme(plot.title=element_text(size=12)) +
          labs(title=NULL, x=NULL, y=NULL) +
          theme(legend.position="none")
        
        eval(parse(text=paste0(
          "p.", type, "[[p.id]] <- p"
        )))
        ggsave(filename=paste0(out.dir, "/", chr, "_", gcb, "_", ijset, "_",
                               p.id, "CsVsCII.pdf"),
               height=10, width=10, unit="in", plot=p); rm(p, p.id)
      } # names(y.v) for loop end
    } # type.v for loop end
    
  }
  
  rm(CSCPCII.MX); gc()
  
} # ct.v for loop end

if(CsCpPlot){
  p.arr <- ggarrange(plotlist=p.lst, nrow=4, ncol=6, legend=NULL)
  ggexport(p.arr, height=40, width=60,
           filename=paste0(out.dir, "/", chr, "_", gcb, "_", ijset, "_cuts", 
                           cuts, "_bins", bins, "_", Cs.form, "CsVsCp.pdf"))
}
if(CsCIIPlot){
  for(type in type.v){
    eval(parse(text=paste0(
      "p.arr <- ggarrange(plotlist=p.", type, ", nrow=4, ncol=6, legend=NULL)"
    )))
    #ggexport(p.arr, height=40, width=60,
    ggexport(p.arr, height=20, width=30,
             filename=paste0(out.dir, "/", chr, "_", gcb, "_", type, "_", ijset, 
                             "_cuts", cuts, "_bins", bins, "_CsVsCII.pdf"))
    rm(p.arr)
  }
}

# rm(list=ls()); gc()
