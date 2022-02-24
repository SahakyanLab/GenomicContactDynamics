################################################################################
# Calculate per class/family/subfamily content in genome
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/17_RepeatAge")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
repeatmx.dir = paste0(wk.dir, "/z_ignore_git/out_hg19Repeats_summary")
out.dir = paste0(wk.dir, "/out_repeatContent")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# REPEAT.MX (base)
load(file=paste0(repeatmx.dir, "/hg19repeats_repName_base.RData"))

rownames(REPEAT.MX) <- NULL
REPEAT.MX <- REPEAT.MX[,c("copyNumber", "repName", "repClass", "repFamily")]

REPEAT.MX$copyNumber <- as.numeric(as.character(REPEAT.MX$copyNumber))
total <- sum(REPEAT.MX$copyNumber, na.rm=F)

col.v <- setdiff(colnames(REPEAT.MX), "copyNumber")

COPYNUM <- list()
for(col in col.v){
  
  agg <- aggregate(x=REPEAT.MX$copyNumber, by=list(REPEAT.MX[[col]]),  
                   FUN=sum, na.rm=F)
  colnames(agg) <- c(col, "copyNumber")
  
  if( sum(agg[,"copyNumber"], na.rm=F) != total ){
    stop(paste0(col, ": Wrong total copyNumber."))
  } else {
    
    agg <- agg[order(agg[,"copyNumber"], decreasing=T),]
    agg$fraction <- agg$copyNumber / total
    
    print( sum(agg$fraction), quote=F)
    
    #if( sum(agg$fraction != 1) ){
    #  stop(paste0(col, ": Fraction does not sum up to 1."))
    #} else {
      
      COPYNUM[[col]] <- agg
      rm(agg)
      
    #}
    
    print(col, quote=F)
    
  }
  
}

COPYNUM$totalcopyNumber <- total

save(COPYNUM, file=paste0(out.dir, "/hg19repeats_copyNumber.RData"))

# rm(list=ls()); gc()
