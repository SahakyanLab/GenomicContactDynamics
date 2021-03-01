################################################################################
# Function to get affix from maskfile
################################################################################
################################################################################
getMaskAffix <- function(mask="ct_ESC_foi_ATF2_desc_TF;ct_ESC_foi_ASH2L_desc_TF",
                         splitChar=";"){
  
  mask <- unlist( strsplit(x=mask, split=splitChar) )
  
  # Isolate affix for naming outputs
  affix <- sapply(X=mask, simplify=TRUE, FUN=function(x){
    strsplit(x=x, split="_desc_", fixed=TRUE)[[1]][1]
  }) 
  affix <- gsub(x=affix, pattern="ct\\_|foi\\_", replacement="")
  affix <- ifelse(length(mask)>1, paste(affix, collapse="_"),affix)
  
  return( paste0("_", affix) )
}
