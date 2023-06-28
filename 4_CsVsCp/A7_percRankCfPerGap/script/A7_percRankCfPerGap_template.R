################################################################################
# Per chr, per cell/tissue, save object with percentile rank of contact based on
# cf at certain gap (excluding contacts with zero cf).
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
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/4_CsVsCp")
melt.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm")
out.dir = paste0(wk.dir, "/out_percRankCfPerGap")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chrarr1.repl" 
ct = "arr2.repl"  
min.gap.bin = 50 # j - i - 1 
bin.len = 40000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrlen.df <- read.delim(file=chrlen.file, stringsAsFactors=F) 
chr.tot.bin <- ceiling(chrlen.df$length.bp[chrlen.df$chromosome == chr] / bin.len)

# Calculate expected number of contacts per chr given min.gap.bin
chr.exp.tot.ij <- ((chr.tot.bin * chr.tot.bin) - chr.tot.bin) / 2

# Load MELT.MX and check if it contains expected number of contacts for chr (all gaps)

load(paste0(melt.dir, "/human_", chr, "_allcontacts.RData"))

melt.tot.ij <- length(MELT.MX$upper.tri$i) + length(MELT.MX$upper.tri.nocontact$i)
if(melt.tot.ij != chr.exp.tot.ij){
  rm(MELT.MX)
  stop(paste0(chr, ": MELT.MX does not contain expected chr total contacts."))
}

# # Get count of contacts not present in any tissue per gap, add these values when
# # calculating expected cf per gap because 0s should be considered. 
# # Expected value should include 0s to reflect how hard it is to form contacts
# # at a gap and to highlight contacts that form despite difficult conditions. 
# # E.g. at a long contact gap in which a contact is very hard to maintain, 
# # imagine a case where there is no contact except for 1. That 1 contact should be 
# # deemed important for forming.
# 
# gaps.no.ij.allct <- MELT.MX$upper.tri.nocontact$j - MELT.MX$upper.tri.nocontact$i - 1 
# rm(MELT.MX)
# gaps.no.ij.allct <- gaps.no.ij.allct[gaps.no.ij.allct >= min.gap.bin]
# gaps.no.ij.allct <- table(gaps.no.ij.allct)

# Get observed HiCNorm cf 

load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
mx <- cbind(obs.val=PERSIST.MX$hits[[ct]], 
            gap=(PERSIST.MX$hits$j - PERSIST.MX$hits$i - 1))
rm(PERSIST.MX)

mx <- cbind(mx, value=NA)

# Get all gaps to consider

gaps.chr <- min.gap.bin:(chr.tot.bin - 1 - 1)

# Calculate observed / expected per gap.val

is.nonZeroCf <- mx[,"obs.val"] > 0

for(gap.val in gaps.chr){
  
  gap.str <- as.character(gap.val) 
    
  is.gap <- mx[,"gap"] == gap.val
  
  # Convert non-zero Cf contacts to percentile rank
  obs.nonZeroij <- mx[is.nonZeroCf & is.gap, "obs.val"]
  
  obs.nonZeroij.rnk <- rank(obs.nonZeroij, ties.method="max")
  obs.nonZeroij.perc.rnk <- obs.nonZeroij.rnk / length(obs.nonZeroij) * 100
  
  mx[is.nonZeroCf & is.gap, "value"] <- obs.nonZeroij.perc.rnk
  
  #
  
  rm(gap.str, is.gap, obs.nonZeroij, obs.nonZeroij.rnk, obs.nonZeroij.perc.rnk)
  message(paste0(chr, " ", gap.val, " gap: done!"))
  
}

zeroCf.len <- sum(mx[,"obs.val"] == 0)
valueNA.len <- sum(is.na(mx[,"value"]))

if(zeroCf.len == valueNA.len){
  save(mx, file=paste0(out.dir, "/", ct, "_", chr, "_percRankCfPerGap_", gcb, ".RData"))
} else {
  stop("Error.")
}

# rm(list=ls()); gc()