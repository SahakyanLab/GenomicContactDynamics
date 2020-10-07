################################################################################
# Calculate confusion matrix-derived metrics
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
csv.dir = paste0(wk.dir, "/out_compare/csv_chr1_whole_maskMidSquare_gap50up_refCp")
out.dir = paste0(wk.dir, "/out_cutoffPlotPerMetric")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1" 
ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
         "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
metric.v = c(subj="SIM.3.2.kmer.5", ref="Cp")
src.id = "whole_maskMidSquare_gap50up" 
out.id = "partRange"

# Closed ranges where most contacts belong; to colour points differently
s.range = 'subj.range = c(-0.0001,0.004)'
r.range = 'ref.range = NULL' #c(-0.02,1)'

confMxMetric.v = c("SPN", "RPN", "TPR", "TNR", "PPV", "NPV",  "FNR", "FPR", "FDR",
                   "FOR", "PT", "TS", "ACC", "BA", "F1", "MCC", "FM", "BM", "MK")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(gridExtra)
library(viridis)
library(ggpubr)
library(RColorBrewer)
source(paste0(wk.dir, "/lib/confusionMxMetric.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
id <- paste(paste(names(metric.v), metric.v, sep="_"), collapse="_")
run.name <- paste(gcb, chr, src.id, id, out.id, sep="_")
print(paste0(run.name, "..."), quote=FALSE)
print(s.range, quote=FALSE)
print(r.range, quote=FALSE)

eval(parse(text=s.range))
eval(parse(text=r.range))

margin.p <- list(geom_point(fill="black"), 
                 labs(x="Cut-off", y="Fraction of contacts"),
                 theme(panel.background=element_rect(colour="gray22", 
                                                     size=1, fill=NA),
                       panel.grid.major=element_blank(),       
                       panel.grid.minor=element_blank(),
                       plot.title=element_text(size=10))
)

p1.lst <- frSUBJ <- frREF <- list()
ct.v.len <- length(ct.v)
for(i in 1:ct.v.len){
  
  ct <- ct.v[i]
  out.name <- paste(gcb, chr, ct, src.id, id, sep="_")
  print(paste0(out.name, "_", out.id, "..."), quote=FALSE)
  
  COMPIJMX <- read.csv(file=paste0(csv.dir, "/", out.name, ".csv"), header=TRUE, 
                       stringsAsFactors=FALSE)
  
  # Percentage of contacts considered during comparison
  percij <- COMPIJMX$final.nonNA.ij/(COMPIJMX$final.NA.ij+COMPIJMX$final.nonNA.ij)
  percij <- format(unique(percij)*100, digits=4) 
  
  # Filter based on subj and ref cut-off ranges
  if( is.null(subj.range ) ){
    subj.range <- c( min(COMPIJMX$c.offsubj), max(COMPIJMX$c.offsubj) )
  }
  if( is.null(ref.range) ){
    ref.range <- c( min(COMPIJMX$c.offref), max(COMPIJMX$c.offref) )
  }
  
  out.range <- (COMPIJMX$c.offsubj<subj.range[1] | COMPIJMX$c.offsubj>subj.range[2]) |
    (COMPIJMX$c.offref<ref.range[1] | COMPIJMX$c.offref>ref.range[2])
  COMPIJMX <- COMPIJMX[out.range==FALSE,]
  rm(out.range)
  
  frSUBJ[[ct]] <- cbind.data.frame(ct=ct, cutoff=COMPIJMX$c.offsubj, 
                                   fr=COMPIJMX$SP/COMPIJMX$final.nonNA.ij)
  frREF[[ct]] <- cbind.data.frame(ct=ct, cutoff=COMPIJMX$c.offref, 
                                  fr=COMPIJMX$RP/COMPIJMX$final.nonNA.ij)
  
  # Plot
  p.lst <- list()
  m.len <- length(confMxMetric.v)
  for(m in 1:m.len){
    
    metric <- confMxMetric.v[m]
    densval <- confMxMetric(CONFMX=COMPIJMX[,c("TP","FP", "TN", "FN", 
                                               "SP", "SN", "RP", "RN"),],
                            metric=metric)
    # Plot, source: https://www.rayshader.com/reference/plot_gg.html
    p.lst[[metric]] <- ggplot(data=cbind(COMPIJMX[,c("c.offsubj", "c.offref")], densval),
                              aes(x=c.offsubj, y=c.offref, z=densval)) + 
      geom_tile(aes(fill=densval)) +
      geom_contour(color="cyan") +
      scale_fill_viridis() + 
      labs(x="Cut-off subj", y="Cut-off ref", fill=metric) +
      
      theme(panel.background=element_rect(colour="gray22", size=1, fill=NA),
            panel.grid.major=element_blank(),       
            panel.grid.minor=element_blank())
    
    if(metric=="MCC"){
      p.top.subj <- ggplot(COMPIJMX, aes(x=c.offsubj, y=SP/final.nonNA.ij)) +
        labs(title=paste0(out.name, "_", percij, "%possibleij_", metric)) + 
        margin.p
      p.right.ref <- ggplot(COMPIJMX, aes(x=c.offref, y=RP/final.nonNA.ij)) +
        margin.p + 
        coord_flip() 
      p.main <- p.lst[[metric]] + guides(fill=FALSE) 
      p1.lst[[ct]] <- arrangeGrob(p.top.subj, p.lst[[metric]], p.main, p.right.ref,
                                 ncol=2, nrow=2, widths=c(3,1), heights=c(1,3))
      
      rm(p.top.subj, p.right.ref, p.main)
    }
    
    p.lst[[metric]] <- p.lst[[metric]] + bgr2
    
    print(paste0(ct, "_", metric, " done!"), quote=FALSE)
    
  }
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=5, ncol=5)
  ggexport(p.arr, height=50, width=50,
           filename=paste0(out.dir, "/", out.name, "_", out.id, ".pdf" ))
  
  rm(p.arr, p.lst, out.name, COMPIJMX)
}

out.name <- run.name
p.arr <- ggarrange(plotlist=p1.lst, nrow=5, ncol=5)
ggexport(p.arr, height=50, width=50, filename=paste0(out.dir, "/", out.name, "_fr.pdf" ))

# Plot fraction of contacts
frSUBJ <- do.call("rbind.data.frame", frSUBJ)
frREF <- do.call("rbind.data.frame", frREF)
coul <- colorRampPalette(rev(brewer.pal(n=11,name="Spectral")))(ct.v.len)
p.lst <- list()
for(x in c("SUBJ", "REF")){
  eval(parse(text=paste0(
    "p.lst[[x]] <- ggplot(data=fr", x, ", aes(x=cutoff, y=fr))"
  )))
  p.lst[[x]] <- p.lst[[x]] +
    geom_point(aes(col=ct), alpha=0.5, size=2.5, shape=1) +
    scale_colour_manual(values=coul) +
    labs(title=paste0(out.name, "_", x), 
         x="Cut-off", y="Fraction of contacts") +
    bgr2
}
p.arr <- ggarrange(plotlist=p.lst, nrow=1, ncol=2)
ggexport(p.arr, height=10, width=20, filename=paste0(out.dir, "/", out.name, "_frCombined.pdf" ))

# rm(list=ls()); gc()

