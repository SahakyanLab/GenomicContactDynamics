################################################################################
# Plot calculated observed over expected Cf vs. Cp.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/4_CsVsCp")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm")
src.dir = paste0(wk.dir, "/out_obsOverExpCfVsCp")
out.dir = paste0(wk.dir, "/out_obsOverExpCfVsCp_plot")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chrs = paste0("chr", c(1:22, "X")) #paste0("chr", c(1:22, "X"))
nCPU = 1
cts = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
             "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "LC", "FC"))
min.gap.bin = 50 # j - i - 1 
gap.rng.closed = c(50, Inf)
bin.len = 40000
Cps = 1:21
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
library(ggplot2)
library(cowplot)
source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(gcb, "_gap.rng.closed", gap.rng.closed[[1]], "To", gap.rng.closed[[2]],
                   "_log2obsOverExpCfVsCp")

chrs.len <- length(chrs)

P.LST <- list()
for(ct in cts){
  
  #### PARALLEL EXECUTION #########
  MX <- foreach(itr=isplitVector(1:chrs.len, chunks=nCPU), 
                .inorder=T, .combine="rbind"
                
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=F, FUN=function(i){
      
      chr <- chrs[[i]]
      load(paste0(src.dir, "/", ct, "_", chr, "_ObsOverExp_", gcb, ".RData"))
      load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      
      if( length(mx[,1]) == length(PERSIST.MX$ntis) ){
        mx <- cbind(mx, Cp=PERSIST.MX$ntis)
      } else {
        stop(paste0(chr, " ", ct, ": Checkpoint 1."))
      }
      
      message(paste0(chr, " ", ct, " done!"))
      return(mx)
      
    })
    
    return(do.call("rbind", chunk))
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  df <- as.data.frame(MX)
  rm(MX)
  
  # Filter 
  df <- df[!is.nan(df$OE),] # obs / exp of NaN means exp.val of 0 so no contact at gap
  df <- df[df$OE != 0,] # Contact not in tissue
  df <- df[ df$gap >= gap.rng.closed[[1]] & df$gap <= gap.rng.closed[[2]], ]
  
  df$Cp <- factor(df$Cp, levels=as.character(Cps))
  
  # Log2 Observed over expected
  df$OE <- log2(df$OE)
  
  P.LST[[ct]] <- ggplot(data=df, aes(x=Cp, y=OE)) +
    stat_boxplot(geom="errorbar", width=0.5) + 
    geom_boxplot(fill="honeydew3", outlier.color=adjustcolor("black", alpha=0.05)) +
    bgr1 + 
    theme(axis.text.x=element_blank(), 
          axis.title.x=element_blank(), axis.title.y=element_blank())
  
  rm(df)
  message(paste0(ct, " done!"))
  
} # cts for loop end

p.arr <- cowplot::plot_grid(plotlist=P.LST, ncol=7, align="hv", byrow=T)
save_plot(filename=paste0(out.dir, "/", out.name, ".png"), plot=p.arr, 
          base_height=3 * 5, base_width=7 * 5)


# rm(list=ls()); gc()