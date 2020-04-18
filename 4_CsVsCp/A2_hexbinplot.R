################################################################################
# Hexbin plot of Cs vs. Cp per cell/tissue.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/13_CsVsCp"
    data.dir = "/Users/ltamon/Database"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/13_CsVsCp"
    data.dir = "/t1-data/user/ltamon/Database"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_hexbin")
### OTHER SETTINGS #############################################################
# Expands warnings
options(warn=1)
ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
         "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "LC", "FC")
gcb = "min05Mb"
#chr.v = paste("chr", c(1:22, "X"), sep="")
chr.v = "chrALL"
nCPU = 1L #~15G
# Scaled Cs values?
scaled = FALSE
approach = "gghexbin" # "hexbin" "gghexbin"
# If approach = "gghexbin"
cuts = 8
plotOnly = TRUE
# Combine plots?
combine = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
library(foreach)
library(doParallel)
library(itertools)
library(viridis)
library(hexbin)
library(ggplot2)
library(ggpubr)
# Loaded only when needed
#library(Hmisc)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/makeHexbinggplot.R"))
### FUNCTION ###################################################################
makeCpVsCsHexbin <- function(
  out.dir = "/dir",
  persist.dir = "/dir",
  ct = "FC",
  gcb = "min2Mb",
  chr.v = paste("chr", c(1:22, "X"), sep=""),
  nCPU = 5L,
  scaled = FALSE,
  plotOnly = FALSE
){
  
  affix <- ifelse(scaled==TRUE, "_scaled", "")
  id <- paste0("chrAll_", gcb, "_", ct)
  
  if(plotOnly==FALSE){
    
    if(nCPU > 1){
      registerDoParallel(cores=nCPU)
      `%op%` <- `%dopar%`
      print(paste0("Running with ", nCPU, " cores."), quote=F)
    } else {
      `%op%` <- `%do%`
    }
    
    chr.v.len <- length(chr.v)
    toExport <- c("affix", "ct", "chr.v", "gcb")
    
    #### PARALLEL EXECUTION #########
    CPCS.MX <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                       .inorder=TRUE, .combine="rbind",
                       .export=toExport, .noexport=ls()[!ls()%in%toExport]
                       
    ) %op% {
      
      chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
        
        chr <- chr.v[i]
        # Load PERSIST.MX (original/scaled)
        load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, affix, ".RData"))
        rownames(PERSIST.MX$hits) <- NULL
        
        log.ct <- PERSIST.MX$hits[,ct]!=0
        
        if( sum(log.ct)>nrow(PERSIST.MX$hits) ){
          stop(paste0(chr, ":Checkpoint 1."))
        }
        
        print(paste0(chr, " done!"), quote=FALSE)
        
        return(
          cbind(Cs=PERSIST.MX$hits[log.ct,ct], 
                Cp=PERSIST.MX$ntis[log.ct])
        )
        
      })
      
      return( do.call("rbind", chunk) ) 
      
    }
    ### END OF PARALLEL EXECUTION ###
    
    save(CPCS.MX, file=paste0(out.dir, "/", id, "_CsVsCp_hexbinplot.RData"))
    
  } else {
    
    load(file=paste0(out.dir, "/", id, "_CsVsCp_hexbinplot.RData"))
    
  }
  
  affix <- gsub(x=affix, pattern="_", replacement="", fixed=TRUE)
 
  # ggplot-based hexbinplot
  if(approach=="gghexbin"){
    library(Hmisc)
    
    id <- paste0(id, "_cuts", cuts)
    
    p <- makeHexbinggplot(xvar=CPCS.MX[,"Cp"], 
                          yvar=CPCS.MX[,"Cs"], 
                          bins=30, 
                          cuts=cuts,
                          xlab=expression( bold("C"["p"]) ),
                          ylab=bquote( bold("C"["s"]~.(affix)) ),
                          title=id,
                          col=viridis(cuts)
    ) 
    
    ggsave(filename=paste0(out.dir, "/", id, "_CsVsCp_hexbinplot.pdf"),
           units="in", width=10, height=10, plot=p$hexplot)
  }
  
  # Hexbinplot way
  if(approach=="hexbinplot"){
    
    pdf(file=paste0(out.dir, "/", id, "_CsVsCp_hexbinplot.pdf"), 
        width=10, height=10)
    print( hexbinplot(Cs~Cp, data=as.data.frame(CPCS.MX), colramp=viridis(12),
               trans=log, inv=exp, aspect=1, mincnt=1,
               xlab=expression( bold("C"["p"]) ), 
               ylab=bquote( bold("log"["10"]~"(C"["s"]~")"~.(affix)) ),
               main=id)
    )
    dev.off()
  }
  
  rm(CPCS.MX, affix, id); gc()
  
  print(paste0(gcb, "_", ct, " done!"), quote=FALSE)
  
  return(p$hexplot)
}
################################################################################
makeCpVsCsHexbin <- cmpfun(makeCpVsCsHexbin, options=list(suppressUndefined=TRUE))
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
p.lst <- list()

for( i in 1:length(ct.v) ){
  
  p.lst[[i]] <- makeCpVsCsHexbin(
    persist.dir=persist.dir,
    out.dir=out.dir,
    ct=ct.v[i],
    gcb=gcb,
    chr.v=chr.v,
    nCPU=nCPU,
    scaled=scaled,
    plotOnly=plotOnly
  )
  
} 

if(combine){
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=3, ncol=7,
                     legend=NULL)
  ggexport(p.arr, width=60, height=24,
           filename=paste0(out.dir, "/", gcb, "_cuts", cuts, "_hexbinplot.pdf" ))
}


# rm(list=ls())
