################################################################################
# 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/binmx/out_bindata_1perc_HiCNorm"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
meanCp.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binWeightedMeanCp")
mutCp.dir = paste0(wk.dir, "/out_mutCalcPerCp/WT_SEQ_rowSum")
out.dir = paste0(wk.dir, "/out_mutCalcVsMeanCp/WT_SEQ_rowSum")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
src.id = "hg38ToHg19" # "Hg19" | "hg38ToHg19"
mut.v = c("All", "A>C", "A>G", "A>T", "C>A", "C>G", "C>T")
calc.v = c("Tmut", "Nmsite", "TmutDIVNmsite", "Nmsitenorm")
hxbincuts.v = c(13,13,3,13)
names(hxbincuts.v) <- calc.v
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(Hmisc)
library(hexbin)
library(viridis)
library(ggpubr)
source(paste0(lib, "/makeHexbinggplot.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
hxbinCp0 <- list()
hxbinCp <- list()
for(mut in mut.v){
  
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=TRUE)
  
  # MUTCP.DF
  load(file=paste0(mutCp.dir, "/", mut.id, "_", src.id, "_mutCalcPerCp.RData"))
  rownames(MUTCP.DF) <- paste0(MUTCP.DF$chr, ".", MUTCP.DF$bins)
  MUTCP.DF$wmeanCp0 <- NA
  MUTCP.DF$wmeanCp <- NA
  chr.v <- unique(MUTCP.DF$chr)
  
  # Add mean Cp of bins in MUTCP.DF
  for(chr in chr.v){
    
    # BINWMEANCP.DF
    load(file=paste0(meanCp.dir, "/", gcb, "_", chr, "_weightedMeanCp.RData"))
    rownames(BINWMEANCP.DF) <- paste0(chr, ".", BINWMEANCP.DF$bin)
    
    if(!all(rownames(BINWMEANCP.DF)%in%rownames(MUTCP.DF))){
      stop(paste0(chr, ": Checkpoint 1."))
    }
    MUTCP.DF[rownames(BINWMEANCP.DF), "wmeanCp0"] <- BINWMEANCP.DF[rownames(BINWMEANCP.DF), "wmeanCp0"]
    MUTCP.DF[rownames(BINWMEANCP.DF), "wmeanCp"] <- BINWMEANCP.DF[rownames(BINWMEANCP.DF), "wmeanCp"]
    rm(BINWMEANCP.DF); gc()
    
  }
  
  for(calc in calc.v){
    
    hxbincuts <- hxbincuts.v[calc]
    hxbinCp0[[paste0(mut, ".", calc)]] <- makeHexbinggplot(xvar=MUTCP.DF[,"wmeanCp0"],
                                                           yvar=MUTCP.DF[,calc], 
                                                           bins=30, 
                                                           cuts=hxbincuts,
                                                           xlab="wmeanCp0",
                                                           ylab=calc,
                                                           title=paste0(mut.id, "_", src.id, "_", calc, "_hxbincuts", hxbincuts),
                                                           col=viridis(hxbincuts))$hexplot
    
    hxbinCp[[paste0(mut, ".", calc)]] <- makeHexbinggplot(xvar=MUTCP.DF[,"wmeanCp"],
                                                          yvar=MUTCP.DF[,calc], 
                                                          bins=30, 
                                                          cuts=hxbincuts,
                                                          xlab="wmeanCp",
                                                          ylab=calc,
                                                          title=paste0(mut.id, "_", src.id, "_", calc, "_hxbincuts", hxbincuts),
                                                          col=viridis(hxbincuts))$hexplot
                                                          
  }
  
  rm(MUTCP.DF, mut.id, hxbincuts); gc()
  print(paste0(mut, " done!"), quote=FALSE)
  
} # mut.v for loop end



calc.v.len <- length(calc.v)
mut.v.len <- length(mut.v)

p.arr <- ggarrange(plotlist=hxbinCp0, nrow=mut.v.len, ncol=calc.v.len, legend=NULL)
ggexport(p.arr, height=mut.v.len*10, width=calc.v.len*10,
         filename=paste0(out.dir, "/", src.id, "_calcXwmeanCp0.pdf"))
rm(hxbinCp0)

p.arr <- ggarrange(plotlist=hxbinCp, nrow=mut.v.len, ncol=calc.v.len, legend=NULL)
ggexport(p.arr, height=mut.v.len*10, width=calc.v.len*10,
         filename=paste0(out.dir, "/", src.id, "_calcXwmeanCp.pdf"))

# rm(list=ls()); gc()

