################################################################################
# Vioplin plots log10(Cs/Cf) vs. Cp
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

data.dir = paste0(home.dir, "/Database")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
lib = paste0(home.dir, "/DPhil/lib")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/4_CsVsCp")
out.dir = paste0(wk.dir, "/out_morePlots")

cscp.id = "HiCNormCs" # rawCs |  HiCNormCs
cscp.dir = paste0(wk.dir, "/out_CsVsCpobj/", cscp.id)
print(paste0("CSCP.MX directory: ", cscp.dir))
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chrALL"
bin.len = 40000
gap.range.bins.closed = c(50,50) # j - i - 1, closed range, NULL to skip filtering
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(vioplot)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(persist.dir, "/chr21_Persist_min2Mb.RData"))
ct.v <- setdiff(colnames(PERSIST.MX$hits), c("i", "j"))
#ct.v <- c("Lu") # REMOVE

out.name <- paste0(gcb, "_", chr, "_", cscp.id, "VsCp_gaprange_", 
                   gap.range.bins.closed[1], "_", gap.range.bins.closed[2])
print(paste0(out.name, "..."), quote=F)

pdf(file=paste0(out.dir, "/", out.name, "_allct_morePlots.pdf"),
    height=10, width=20)

p.lst <- list()
par(mfrow=c(1,2))
for(ct in ct.v){
  
  load(paste0(cscp.dir, "/", chr, "_", gcb, "_", bin.len, "bpbin_", ct, "_CsVsCp.RData"))
  
  if( !is.null(gap.range.bins.closed) ){
    CPCS.MX <- CPCS.MX[ CPCS.MX[,"gap.jminusiminus1.bin"] >= gap.range.bins.closed[1] & 
                        CPCS.MX[,"gap.jminusiminus1.bin"] <= gap.range.bins.closed[2], ]
    gap.range.bins.derived <- range(CPCS.MX[,"gap.jminusiminus1.bin"])
  } else {
    print(paste0(ct, ": No gap filtering."))
  }
  
  CPCS.MX <- as.data.frame(CPCS.MX)
  CPCS.MX$Cp <- factor(as.character(CPCS.MX$Cp), levels=as.character(1:21))
  
  minval.allCp <- format(min(CPCS.MX$Cs), scientific=3, digits=4)
    
  plot.title <- paste0(out.name, "_", ct, "_minimumCsAllCp", minval.allCp,
                       "_Derivedgaprange_", gap.range.bins.derived[1], "_", gap.range.bins.derived[2])
  
  vioplot(formula=Cs~Cp, data=CPCS.MX, las = 1, plotCentre="line",
          xlab="Cp", ylab="Cf", main=paste0(plot.title, "_horizontalLineIsMedian"), col="#FDC776", cex.main=0.8)
  
  boxplot(formula=Cs~Cp, data=CPCS.MX, outline=F, 
          xlab="Cp", ylab="Cf", main=paste0(plot.title, "_outline=F"), col="#FDC776", cex.main=0.8)

  print(paste0(ct, " done!"))
  
  rm(CPCS.MX)
  gc()
  
}

dev.off()

# rm(list=ls()); gc()