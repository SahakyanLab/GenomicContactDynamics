################################################################################
# Visualise distribution of counts of rates per bin to determine the appropriate
# minimum count a bin should have to be considered. Counts made sure to be
# identical with countPerBin column in FEATURE.BIN.MX.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/24_Recombination")
src.dir = paste0(wk.dir, "/z_ignore_git/out_mapToHiCcontactPersistBins")
out.dir = paste0(wk.dir, "/out_countValuesPerBin")
### OTHER SETTINGS #############################################################
chrs = paste0("chr", c(21:22))
out.id = paste0(chrs[[1]], "To", tail(chrs, n=1))
src.id = "min2Mb_recomRates_2011_01_phaseII_B37_Myers"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(Rmisc)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(out.id, "_", src.id)

# Get count of rates of bins from all chr
COUNTBIN <- sapply(X=chrs, simplify=F, FUN=function(chr){
  
  # FEATURE.BIN.MX only contains bins forming a long-range contact
  load(paste0(src.dir, "/", chr, "_", src.id, ".RData"))
  rates <- FEATURE.BIN.MX[,"V5"]
  rates <- rates[!is.na(rates)] 
  rates <- strsplit(x=rates, split=";")
  
  countPerBins <- FEATURE.BIN.MX[,"countPerBin"]
  countPerBins <- countPerBins[countPerBins > 0 & !is.na(countPerBins)]
  countPerBins.derived <- lengths(rates)
  
  if( !identical(as.numeric(countPerBins), as.numeric(countPerBins.derived)) ){
    rm(rates)
    stop( paste0(chr, ": FEATURE.BIN.MX countPerBin not identical to countPerBins.derived
                 derived from FEATURE.BIN.MX rates column." ) )
  } else {
    return( lengths(rates) )
  }
  
})
COUNTBIN <- do.call("c", COUNTBIN)
COUNTBIN <- data.frame(count=COUNTBIN)
rownames(COUNTBIN) <- NULL

# Summary statistics

SUMM.STAT <- summarySE(COUNTBIN, measurevar="count", .drop=F)
save(SUMM.STAT, file=paste0(out.dir, "/", out.name, "_CountOfRatesPerNonNABin_summstat.RData"))

# Plot
p <- ggplot(data=COUNTBIN, aes(x=count)) +
  geom_density(col="#55bde6", fill="#55bde6") +
  labs(title=paste0(out.name, "_CountOfRatesPerNonNABin")) + 
  bgr1

ggsave(filename=paste0(out.dir, "/", out.name, "_CountOfRatesPerNonNABin_density.pdf"),
       height=10, width=10, plot=p)

# rm(list=ls()); gc()