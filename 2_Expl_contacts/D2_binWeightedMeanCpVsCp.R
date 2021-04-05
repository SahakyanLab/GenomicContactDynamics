################################################################################
# Plot weighted mean Cp vs. Cp of bins because a bin with a high Cp may have a
# relatively low weighted mean Cp if it participates in many low Cp contacts
# and very few high Cp contacts. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/2_Expl_contacts"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_Expl_contacts"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
meanCp.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binWeightedMeanCp")
out.dir = paste0(wk.dir, "/out_binWeightedMeanCpVsCp")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X")) 
Cp.v = 1:21
bin.len = 4e4
hxbincuts = 8
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(Hmisc)
library(hexbin)
library(viridis)
library(ggpubr)
source(paste0(lib, "/makeHexbinggplot.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
df <- list()
for(chr in chr.v){
  
  # BIN.MX
  load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
  BIN.MX <- BIN.MX[,!colnames(BIN.MX)%in%c("start","end")]
  
  # BINWMEANCP.DF
  load(file=paste0(meanCp.dir, "/", gcb, "_", chr, "_weightedMeanCp.RData"))
  
  if( !identical(rownames(BIN.MX), as.character(BINWMEANCP.DF$bin)) ){
    stop(paste0(chr, ": Checkpoint 1."))
  }
  
  # Assign bin to Cp (involved in at least 1 contact for that Cp in any tissue)
  for(Cp in Cp.v){
    
    Cp.ind <- grep(x=colnames(BIN.MX), pattern=paste0("s_Cp_", Cp, "_ct"))
    if(length(Cp.ind)!=length(Cp.v)){ stop(paste0(chr, ": Checkpoint 2.")) }
    BINWMEANCP.DF[[as.character(Cp)]] <- as.numeric(rowSums(x=BIN.MX[,Cp.ind])>0)
    rm(Cp.ind)
    
  }
  
  df[[chr]] <- cbind.data.frame(chr=chr, BINWMEANCP.DF, stringsAsFactors=FALSE)
  
  rm(BIN.MX, BINWMEANCP.DF); gc()
  
  print(paste0(chr, " done!"), quote=FALSE)
  
} # chr.v for loop end

df <- do.call("rbind.data.frame", df)
rownames(df) <- NULL

df <- reshape2::melt(data=df, id.vars=c("chr", "bin", "wmeanCp0", "wmeanCp"))
# Choose only bins with Cp (only bins forming long-range contacts)
df <- df[df$value==1,]

#  Boxplot

pdf(file=paste0(out.dir, "/", gcb, "_meanCpVsCp.pdf"), height=10, width=10)
par(mfrow=c(2,2))

hx.lst <- list()
for( wm in c("wmeanCp0", "wmeanCp") ){
  
  eval(parse(text=paste0(
    'boxplot(', wm, '~variable, outline=TRUE, data=df, xlab="Cp", ylab="', wm, '", 
            boxwex=0.6, cex.axis=1.2, col="#FDC776", cex.main=1, 
            main=paste0(gcb, "_', wm, 'VsCp"))'
  )))
  
  eval(parse(text=paste0(
    'boxplot(', wm, '~variable, outline=FALSE, data=df, xlab="Cp", ylab="', wm, '", 
    boxwex=0.6, cex.axis=1.2, col="#FDC776", cex.main=1, 
    main=paste0(gcb, "_', wm, 'VsCp"))'
  )))
  
  hx.lst[[wm]] <- makeHexbinggplot(xvar=df$variable,
                                   yvar=df[[wm]], 
                                   bins=30, 
                                   cuts=hxbincuts,
                                   xlab="Cp",
                                   ylab=wm,
                                   title=paste0(gcb, "_", wm, "VsCp_hxbincuts", hxbincuts),
                                   col=viridis(hxbincuts))$hexplot

}

dev.off()

p.arr <- ggarrange(plotlist=hx.lst, nrow=1, ncol=2, legend=NULL)
ggexport(p.arr, height=10, width=20,
         filename=paste0(out.dir, "/", gcb, "_meanCpVsCp_hxbincuts", hxbincuts, ".pdf"))


# rm(list=ls()); gc()




