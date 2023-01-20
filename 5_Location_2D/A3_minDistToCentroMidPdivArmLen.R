################################################################################
# Combining all chromosome data, plot minimum distance of contact bin midpoint
# to centromere midpoint divided by length of arm, where bin is located. If
# value is positive (because min distance is positive) the bin is located downstream
# of centromere midpoint, if the value is negative (because min distance is negative) 
# the bin is located upstream of centromere midpoint. Because the minimum distance
# values were calculated as bin midpoint minus centromere midpoint using 1-based
# coordinates, the same is done for arm lengths i.e. centromere midpoints minus
# chromosome start or end in 1-based CS.
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
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/5_Location_2D")
src.dir = paste0(wk.dir, "/out_distToRefPoints")
out.dir  = paste0(wk.dir, "/out_minDistToCentroMidPdivArmLen")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt") 
# 1-based, constant centromere length of 3Mb
centro.file = paste0(data.dir, "/ucsc_tables/hsa_centromere/ct_hg19_foi_centromereonly_desc_DNA")
# Should match ref.points order in src data (IJDISTTOREF.MX)
ref.points.name = c("chr.start.bp", "chr.end.bp", "centromere.midP.bp") 
### OTHER SETTINGS #############################################################
Cp.v = 1:21
gcb = "min2Mb"
chrs = paste0("chr", c(1:22, "X"))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(RColorBrewer)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Centromere midpoints 
centro.df <- read.delim(file=centro.file, header=F, sep="\t", stringsAsFactors=F)
centro.len <- unique( centro.df$V3 - centro.df$V2 + 1 ) # = 3e+06 bp
if( any(centro.len <= 0) ){
  stop("Negative or zero centromere length.")
}
centro.points <- centro.df$V2 + ((centro.len - 1) / 2) # R rounding off

# Chromosome lengths
chrlen.df <- read.delim(file=chrlen.file, header=T, sep="\t", stringsAsFactors=F)
chrlen.df <- chrlen.df[chrlen.df$chromosome != "chrMT",]

# Chromosome arm lengths
if( identical(centro.df$V1, chrlen.df$chromosome) ){
  
  arm.df <- matrix(data=NA, nrow=length(chrlen.df[,1]), ncol=2, dimnames=list(chrlen.df$chromosome, NULL))
  arm.df[,1] <- centro.points - 1
  arm.df[,2] <- chrlen.df$length.bp - centro.points
  arm.df <- as.data.frame(arm.df)

} else {
  stop("Order of chromosomes between centro.df and chrlen.df not the same")
}

# Plot 

ref.points.id <- paste(ref.points.name, collapse="_")
DF <- sapply(X=chrs, simplify=F, FUN=function(chr){
  
  load(paste0(src.dir, "/", gcb, "_", chr, "_ijdisttoref_", ref.points.id, ".RData"))
  df <- as.data.frame(IJDISTTOREF.MX, stringsAsFactors=F)
  df$chr.start.bp <- df$chr.end.bp <- NULL
  rm(IJDISTTOREF.MX)
  
  # Determine which arm the bin is located
  df$arm <- NA
  df$arm[df$centromere.midP.bp < 0] <- 1
  df$arm[df$centromere.midP.bp > 0] <- 2
  binInCentroMidP.len <- sum(is.na(df$arm))
  if(binInCentroMidP.len > 0){
    warning(paste0(chr, ": ", binInCentroMidP.len, " row/bin midpoint/s distance to centromere midpoint equal 0.s"))
    df$arm[is.na(df$arm)] <- 1
  }
  
  # Get fraction value of distance from centromere midpoint relative to length of arm location
  arm.lens <- as.numeric(arm.df[chr,])
  armlen.divisor <- arm.lens[as.numeric(df$arm)]
  df$minDcentmidDivArm <- df$centromere.midP.bp / armlen.divisor
  df$centromere.midP.bp <- df$arm <- NULL 
  
  message(paste0(chr, " done!"))
  return(df)
  
})

DF <- do.call("rbind", DF)
rownames(DF) <- NULL

DF$Cp <- factor(as.character(DF$Cp), levels=as.character(Cp.v))
coul <- colorRampPalette( rev( brewer.pal(11, "Spectral") ) )( length(levels(DF$Cp)))

p <- ggplot(data=DF, aes(x=minDcentmidDivArm, group=Cp)) +
  geom_density( position="identity", aes(colour=Cp) ) + 
  scale_colour_manual(values=coul) + 
  labs(title=paste0(gcb, "_chrALL_minDcentmidDivArmLength_asindistToRefPoints_alldistareCoordA-CoordB_1based
                    _mindistwith0dividedbyarmlen1but0Div1isStill0.Didthistoseeifnohighdensityat0."), 
       x="minDcentmidDivArmLength", colour="Cp") +
  guides(colour=guide_legend(ncol=1)) +
  bgr2 

ggsave(filename=paste0(out.dir, "/", gcb, "_minDcentmidDivArm_", ref.points.id, "_density.pdf"),
       units="in", width=10, height=10, plot=p)

# rm(list=ls()); gc()