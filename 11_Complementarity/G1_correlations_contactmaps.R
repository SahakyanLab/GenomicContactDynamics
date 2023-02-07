################################################################################
# Correlations between Cp, Cf and CII per chromosome, per tissue per gap subranges
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity")
melt.id = "HiCNorm_primary_cohort"
melt.dir = paste0(home.dir, "/Database/GSE87112/combined_contacts/", melt.id)
cii.dir = paste0(wk.dir, "/out_constraints_GfreeSingleNorm/merged_final")
out.dir = paste0(wk.dir, "/out_correlations_contactmaps")
### OTHER SETTINGS #############################################################
cii.id = "kmer_min2Mb"
chr = "chr17"
ct = "FC"
bin.len = 40000
gap.res.bp = 5e6
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(data.table)
library(viridis)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/doCorTest.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(chr, "_", ct, "_gap.res.bp", (gap.res.bp / 1e6), "Mb_bin.len", 
                   bin.len, "bp_", melt.id, "_NonZeroCfOnly", cii.id)

# Load MELT.MX
load(paste0(melt.dir, "/human_", chr, "_allcontacts.RData"))
# Load CII.MX
load(paste0(cii.dir, "/", chr, "_", cii.id, ".RData"))

# Append Cf data to CII.MX

wrong.TF <- c(F,F)
wrong.TF[1] <- !identical( as.integer(CII.MX[rownames(MELT.MX$upper.tri),"i"]),
                           as.integer(MELT.MX$upper.tri[,"i"]) )

wrong.TF[2] <- !identical( as.integer(CII.MX[rownames(MELT.MX$upper.tri),"j"]),
                           as.integer(MELT.MX$upper.tri[,"j"]) )

if( any(wrong.TF) ){
  rm(CII.MX)
  stop(paste0(chr, ": Error when appending Cf data to CII.MX"))
}

MX <- cbind(CII.MX, Cf=NA)
MX[rownames(MELT.MX$upper.tri),"Cf"] <- MELT.MX$upper.tri[,ct]

rm(CII.MX, MELT.MX)

# Identify gap subranges

gaps <- unname(MX[,"j"] - MX[,"i"] - 1)
gap.limits <- range(gaps)
gap.res.bin <- gap.res.bp / bin.len
gap.subrng.lowerb <- seq(gap.limits[1], gap.limits[2], gap.res.bin)
gap.subrng.mx <- cbind(gap.subrng.lowerb,
                       c(gap.subrng.lowerb[-1], gap.limits[2])
                       )
# Remove 0 length subranges when lower bound maximum is equal to gap maximum
gap.subrng.mx <- gap.subrng.mx[ (gap.subrng.mx[,2] - gap.subrng.mx[,1]) > 0, ]
gap.subrng.len <- length(gap.subrng.mx[,1])

if( gap.subrng.len != ceiling(gap.limits[2] / gap.res.bin) ){
  rm(gap.subrng.mx)
  stop(paste0(chr, ": gap.subrng.len not equal to expected number of gap subranges given bin.len and maximum gap or gap.limits[2]."))
}

# Correlations per gap subrange
P.LST <- COR.OUT <- list()
for(i in 1:gap.subrng.len){
  
  is.gap.closed <- (gaps >= gap.subrng.mx[i,1]) & (gaps <= gap.subrng.mx[i,2])
  is.inclij <- MX[,"Cf"] > 0 & !is.na(MX[,"Cf"]) & !is.na(MX[,"C||"]) & is.gap.closed
  
  COR.OUT[[i]] <- doCorTest(xval=as.numeric(MX[is.inclij,"Cf"]), yval=as.numeric(MX[is.inclij,"C||"]),
                            out.dir=NULL, out.name=NULL)
 
  # Plot points
  
  df <- as.data.frame(MX[is.inclij,c("Cf", "C||")])
  setnames(df, old="C||", new="CII")
  
  plot.title <- paste0(paste0("GapSubrange", i), 
                       "\nPearCoef_", COR.OUT[[i]]$pear$estimate, "_pval", COR.OUT[[i]]$pear$p.value,
                       "\nSpearCoef_", COR.OUT[[i]]$spea$estimate, "_pval", COR.OUT[[i]]$spea$p.value)
  
  P.LST[[i]] <- ggplot(data=df, aes(x=Cf, y=CII)) + 
    geom_hex() + 
    stat_binhex(aes(label=..count..), geom="text", colour="white", size=1) + 
    scale_fill_viridis() +
    labs(title=plot.title) + 
    bgr1 +
    theme(plot.title=element_text(size=20))
  
  rm(is.gap.closed, is.inclij, df)
  
}

save(COR.OUT, file=paste0(out.name, "_OrderedByGapSubrange_cortest.RData"))

p.arr <- ggarrange(plotlist=P.LST, nrow=2, ncol=3, legend=NULL)
ggexport(p.arr, height=20, width=30, filename=paste0(out.dir, "/", out.name, "_CfvsCII_hex.pdf"))


# rm(list=ls()); gc()