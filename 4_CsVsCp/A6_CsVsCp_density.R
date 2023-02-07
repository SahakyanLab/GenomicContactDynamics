################################################################################
# Generate a across-dataset normalised Cf vs. Cp plot combining data from all
# cell lines/tissues and chromosomes for a given gap range. All contacts with 
# Cp data are considered. Two version of the plot is made one plots all Cf values
# of a contact with Cp data meaning that it plots the 20 zero Cf of a Cp=1 contact.
# The other version only plots the Cf > 0 values of contacts with Cp data.
# If we are interested in answering whether persistent/dynamic contacts tend to
# have higher/lower Cf value in a tissue, it would be better to exclude the 0s
# because the dynamic contacts will have a lot of 0 Cf pulling down the 
# distribution. Currently A5_morePlots.R considers the 0 Cf values. Both box/violin
# and density plots are made. For the density, zero Cf values are translated by the
# minimum non-zero Cf cause density plots log10 transformed Cf and log10 of 0
# is -Inf.
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/4_CsVsCp")
out.dir = paste0(wk.dir, "/out_CsVsCp_density")
melt.id = "HiCNorm_QQ_primary_cohort"
melt.dir = paste0(home.dir, "/Database/GSE87112/combined_contacts/", melt.id)
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
### OTHER SETTINGS #############################################################
persist.id = "Persist_min2Mb"

Cps = 1:21
chrs ="chr21" 
#cts = NULL # If NULL, use all datasets

regenerateData = T

# Affects plot only, source data for plot contacts all upper.tri contacts from
# MELT.MX$upper.tri
gap.range.bins.closed = c(50,50)
bin.len = 40000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(RColorBrewer)
library(reshape2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(persist.id, "_", melt.id, "_CfVsCp")

if(regenerateData){

  for(chr in chrs){
    
    # Load MELT.MX
    load(paste0(melt.dir, "/human_", chr, "_allcontacts.RData"))
    MELT.MX$rest <- MELT.MX$upper.tri.nocontact <- NULL
    
    # Confirm that Cf >= 0
    ct.nmes <- setdiff(colnames(MELT.MX$upper.tri), c("i", "j"))
    if( any(MELT.MX$upper.tri[,ct.nmes] < 0) ){
      rm(MELT.MX)
      stop(paste0(chr, ": Negative Cf value."))
    }
    
    # Load PERSIST.MX
    load(paste0(persist.dir, "/", chr, "_", persist.id, ".RData"))
    
    # Append Cf data with Cp data
    MELT.MX$upper.tri$Cp <- NA 
    MELT.MX$upper.tri[rownames(PERSIST.MX$hits), "Cp"] <- PERSIST.MX$ntis
    
    df <- melt(MELT.MX$upper.tri, id=c("i", "j", "Cp"))
    df <- df[!is.na(df$value),]
    
    save(df, file=paste0(out.dir, "/", chr, "_", out.id, ".RData"))
    
    rm(PERSIST.MX, MELT.MX)
    
  }
  
}

# Combine chr data for plotting, do gap filtering
DF <- sapply(chrs, simplify=F, FUN=function(chr){
  
  load(paste0(out.dir, "/", chr, "_", out.id, ".RData"))
  df$gap.jminusiminus1.bin <- df$j - df$i - 1

  # Contact gap filtering
  
  if( !is.null(gap.range.bins.closed) ){
    df <- df[ df[,"gap.jminusiminus1.bin"] >= gap.range.bins.closed[1] & 
              df[,"gap.jminusiminus1.bin"] <= gap.range.bins.closed[2], ]
  } else {
    print(paste0(ct, ": No gap filtering."))
  }
  
})
DF <- do.call("rbind", DF)
rownames(DF) <- NULL

## Plot

gap.range.bins.derived <- range(DF$gap.jminusiminus1.bin)
DF$i <- DF$j <- DF$gap.jminusiminus1.bin <- NULL
out.name <- paste0(out.id, "VsCp_gaprange_", gap.range.bins.closed[1], "_", gap.range.bins.closed[2])
plot.id <- paste0(out.name, "_Derivedgaprange_", gap.range.bins.derived[1], "_", gap.range.bins.derived[2])

DF$Cp <- factor(as.character(DF$Cp), levels=as.character(sort(Cps)))
coul <- colorRampPalette( rev(brewer.pal(11, "Spectral")) )(length(levels(DF$Cp)))

# Plot with Cf >= 0

for( toplot.id in c("withZeroCf", "NonZeroCf") ){ # Order of toplot.ids vector should be fixed
  
  out.name.fin <- paste0(out.name, "_", toplot.id)
  
  if(toplot.id == "NonZeroCf"){
    DF <- DF[DF$value > 0, ]
    message("Removed zero Cf.")
  }
  
  # Add p-value calculation
  
  Cf.rng.id <- paste(range(DF$value), collapse="-")
  plot.title <- paste0(plot.id, "_Cfrange_", Cf.rng.id)
  
  p <- ggplot(data=DF, aes(x=Cp, y=value)) +
    geom_violin(width=1, fill="#55bde6", scale="width", col="black", trim=T, position=position_dodge(0.5)) +
    geom_boxplot(lwd=0.6, width=0.3, fill="#55bde6") + 
    scale_x_discrete(limits=levels(DF$Cp)) + 
    labs(x="Cp", y="Cf", title=plot.title) + 
    bgr1 + 
    theme(plot.title=element_text(size=7)) 
  
  ggsave(filename=paste0(out.dir, "/", out.name.fin, "_violin.pdf"),
         width=10, height=10, plot=p)
  
  if(toplot.id == "withZeroCf"){
    translate.min <- min(DF$value[DF$value > 0])
    message("Translate Cf values by minimum non zero Cf for density.")
  } else {
    translate.min <- 0
  }
  
  d <- ggplot(data=DF, aes(log10(value + translate.min))) + 
    geom_density(aes(colour=Cp)) +
    labs(x="log10(Cf)", y="Density", title=plot.title) +
    scale_colour_manual(values=coul, limits=levels(DF$Cp)) +
    guides(colour=guide_legend(ncol=1)) + 
    bgr1 +
    theme(plot.title=element_text(size=7)) 
  
  ggsave(filename=paste0(out.dir, "/", out.name.fin, "_density.pdf"),
         width=10, height=10, plot=d)
  
  rm(p, d)
  message(paste0(toplot.id, " done!"))
  
}
 

# rm(list=ls()); gc()