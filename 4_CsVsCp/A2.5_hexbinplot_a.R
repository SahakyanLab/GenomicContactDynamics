################################################################################
# Simplified version for making hexbin plot of log10 Cf vs. Cp per cell/tissue.
# Cf subtracted by min value per cell type then added by 1.001 then 
# log10 transformed. This was done to have transformed Cf > 0. Counts in legend 
# are in millions. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = paste0("/Users/ltamon")
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/4_CsVsCp")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = paste0("/project/sahakyanlab/ltamon")
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/4_CsVsCp")
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
src.dir = paste0(wk.dir, "/out_CsVsCpobj/HiCNormCs")
out.dir = paste0(wk.dir, "/out_hexbinplot_a")
### OTHER SETTINGS #############################################################
cts = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
             "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "LC", "FC"))
src.id = "chrALL_min2Mb_40000bpbin"
nCPU = 2L 
Cps = 1:21
type.plot = "box" # "hexbin"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
library(foreach)
library(doParallel)
library(itertools)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cts.len <- length(cts)

#### PARALLEL EXECUTION #########
P.LST <- foreach(itr=isplitVector(1:cts.len, chunks=nCPU), .inorder=T, .combine="c"
                   
) %op% {
  
  chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
    
    ct <- cts[[i]]
    
    src.nme <- paste0(src.id, "_", ct)
    load(file=paste0(src.dir, "/", src.nme, "_CsVsCp.RData"))
    if( any(CPCS.MX[,"Cs"] <= 0) ){
      rm(CPCS.MX)
      stop(paste0(src.nme, ": Zero/negative Cs.")) 
    }
    
    df <- as.data.frame(CPCS.MX[,c("Cs", "Cp")])
    rm(CPCS.MX)
    df$Cp <- factor(as.character(df$Cp), levels=as.character(Cps))
    
    # Transform Cs
    df$Cs <- log10( df$Cs - min(df$Cs) + 1.001 )
    
    # Base ggplot
    
    p0 <- ggplot(data=df, aes(x=Cp, y=Cs)) + # REMOVE
      labs(title=src.nme) + 
      bgr1 +
      theme(plot.title=element_blank(), axis.text.x=element_blank(), 
            axis.title.x=element_blank(), axis.title.y=element_blank(), 
            legend.title=element_blank())  
      
    if(type.plot == "hexbin"){
      
      p <- p0 + 
        geom_hex() + 
        stat_binhex(aes(label=after_stat(count)), geom="text", colour="black", size=1) + 
        scale_y_continuous(labels=function(y) sprintf("%.2f", y)) +
        scale_fill_distiller(palette="Spectral", labels=function(fll) sprintf("%.1f", (fll / 1e6) )) 
        #scale_fill_viridis() +
        #scale_fill_gradientn(colors=brewer.pal(n=11,"Spectral"), breaks=seq(1,1000,100))
        
      
    } else if(type.plot == "box"){
      
      p <- p0 +
        stat_boxplot(geom="errorbar", width=0.5) + 
        geom_boxplot(fill="honeydew3", outlier.color=adjustcolor("black", alpha=0.1)) +
        scale_y_continuous(labels=function(y) sprintf("%.2f", y))
      
    } else {
      stop("Invalid type.plot argument")
    }
    
    return(p)
    
  })
  
  return(chunk)
  
}
### END OF PARALLEL EXECUTION ###

p.arr <- cowplot::plot_grid(plotlist=P.LST, ncol=7, align="hv", byrow=T)
save_plot(filename=paste0(out.dir, "/", src.id, "_CsVsCp_", type.plot, "plot.png"), plot=p.arr, 
          base_height=3 * 5, base_width=7 * 5)

# rm(list=ls()); gc()
