################################################################################
# Generate PAIRCOR.MX version setting NA to non-significant correlations. 
# Adjust p-values per chromosome then set correlation of pairs > 0.05 Benjamini-
# Hochberg adjusted to NA
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/23_Hub_Gene_Expression/out_coexpression_pairCor")
src.dir = paste0(wk.dir, "/orig_with_pval_column")
out.dir = wk.dir
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
fles <- list.files(path=src.dir, pattern=".RData")

pairCount.mx <- matrix(data=NA, nrow=length(fles), ncol=3, 
                       dimnames=list(fles, c("orig", "noPvalCorRemoved", "sigOnly")))

for(fle in fles){
  
  load(paste0(src.dir, "/", fle))
  
  pairCount.mx[fle, "orig"] <- sum(!is.na(PAIRCOR.MX[,"cor"]))
    
  # Select non-NA pvalues then adjust
  
  pval <- PAIRCOR.MX[,"pval"]
  is.nonNA.pval <- !is.na(pval)
  adj.noNA.pval <- p.adjust(pval[is.nonNA.pval], method="BH")
  
  #plot(density(PAIRCOR.MX[,"pval"], na.rm=T))
  #lines(density(adj.pval, na.rm=T), col="red")
  #lines(density( adj.noNA.pval, na.rm=T), col="green")
  
  adj.pval <- pval
  adj.pval[ is.nonNA.pval ] <- adj.noNA.pval
  
  # Set to NA cases that cor is finite and pval is not 
  
  is.noPvalCor <- is.finite(PAIRCOR.MX[,"cor"]) & !is.finite(pval)
  #if( sum(is.check) > 0 ){
    #rm(PAIRCOR.MX)
    #stop(paste0(fle, ": Checkpoint 1"))
  #}
  PAIRCOR.MX[is.noPvalCor,"cor"] <- NA
  
  pairCount.mx[fle, "noPvalCorRemoved"] <- sum(!is.na(PAIRCOR.MX[,"cor"]))
  
  # Set insignificant correlations to NA
  
  is.setNA.pair <- is.finite(PAIRCOR.MX[,"cor"]) & is.finite(adj.pval) & (adj.pval > 0.05)
  PAIRCOR.MX[is.setNA.pair,"cor"] <- 0
  
  #
  
  pairCount.mx[fle, "sigOnly"] <- sum(!is.na(PAIRCOR.MX[,"cor"]))
  
  save(PAIRCOR.MX, file=paste0(out.dir, "/", fle))
  rm(PAIRCOR.MX)
  
  message(paste0(fle, " done!"))
  
}

pairCount.mx <- cbind(pairCount.mx, 
                      dropoutPerc=(1 - (pairCount.mx[,"sigOnly"] / pairCount.mx[,"orig"])) * 100)
write.csv(pairCount.mx, file=paste0(out.dir, "/pairCountmx.csv"), row.names=T)

# rm(list=ls()); gc()