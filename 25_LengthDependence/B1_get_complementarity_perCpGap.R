################################################################################
# Per gcb, chr, gap and complementarity method, get values for variable and
# persistent contacts
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar = TRUE) 
options(warn = 1) 

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = file.path(home.dir, "git/GenomicContactDynamics/lib/")
wk.dir = file.path(home.dir, "git/GenomicContactDynamics/25_LengthDependence")
data.dir = file.path(home.dir, "data/gcd/Database")
persist.dir = file.path(data.dir, "HiC_features_GSE87112_RAWpc")
#compl.dir = file.path(home.dir, "data/gcd/11_Complementarity/out_constraints_GfreeSingleNorm/merged_final")
compl.dir = file.path(home.dir, "data/gcd/11_Complementarity/out_constraints/merged_final")
out.dir = file.path(wk.dir, "out_get_complementarity_perCpGap")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chrs = paste0("chr", c("X", 22:1)) #c("chr21", "chr22") #paste0("chr", c("X", 1:22))
compl.type = "align" # kmer | align
var.Cp = 1:3
per.Cp = 19:21
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(tidyr)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for (chr in chrs) {
  
  load(file.path(compl.dir, paste0(chr, "_", compl.type, "_", gcb, ".RData")))
  load(file.path(persist.dir, paste0(chr, "_Persist_", gcb, ".RData")))
  
  # Get only contacts with Cp score
  CII.MX <- CII.MX[rownames(PERSIST.MX$hits),]
  if (!identical(as.integer(PERSIST.MX$ntis), as.integer(CII.MX[,"Cp"]))) {
    stop(gcb, " ", chr, ": Cp in PERSIST.MX and CII.MX not matching")
  }
  rm(PERSIST.MX)
  
  CII.MX <- CII.MX[CII.MX[,"Cp"] %in% c(var.Cp, per.Cp), ]
  
  # Complementarity values per Cp (persistent or variable) per gap
  
  use_cols <- intersect(c("Cp", "C||", "Gfree"), colnames(CII.MX))
  gap_vals <- CII.MX[,"j"] - CII.MX[,"i"] - 1
  gap_factor <- factor(as.character(gap_vals), levels = sort(unique(gap_vals)))
  compl_per_cpgap_lst <- split(as.data.frame(CII.MX[,use_cols]), gap_factor)
  
  for (gap_char in names(compl_per_cpgap_lst)) {
    
    Cp_group_factor <- factor(ifelse(compl_per_cpgap_lst[[gap_char]][,"Cp"] %in% per.Cp, "per.Cp", "var.Cp"),
                              levels = c("var.Cp", "per.Cp"))
    saveRDS(split(compl_per_cpgap_lst[[gap_char]][,"C||"], Cp_group_factor),
            file.path(out.dir, paste0(gcb, "_gap", gap_char, "_", chr, "_", compl.type, "_per_cpgroup.rds")))
    
    if (compl.type == "kmer") {
      
      saveRDS(split(compl_per_cpgap_lst[[gap_char]][,"Gfree"], Cp_group_factor),
              file.path(out.dir, paste0(gcb, "_gap", gap_char, "_", chr, "_Gfree_per_cpgroup.rds")))
      
    }
            
  }
  
  rm(CII.MX, gap_vals, gap_factor, Cp_group_factor, compl_per_cpgap_lst)
  message(gcb, " ", chr, " done!")
  
}

# rm(list=ls()); gc()
