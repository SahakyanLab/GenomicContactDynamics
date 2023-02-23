################################################################################
# Transform bioav_rt_withna.RData to be amenable for subsequent code. 
# rowMeans(rt.mx, na.rm=T) gives NaN when taking mean with all NAs (for median, 
# it returns NA). This NaNs are converted to NAs. Note that unlike old form of
# dataset, the final data does not contain all bins of a chr. This is dealth with
# in A1_plotdata.R.
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
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/20_ReplicationTiming")
out.dir = paste0(wk.dir, "/out_transform_rtdata")
reptimePath = paste0(data.dir, "/replication_timing/TARG/rsc/rt_dataset_process/new/bioav_rt_withna.RData")
### OTHER SETTINGS #############################################################
min.cellline = 59 # see /Users/ltamon/Database/replication_timing/A3_explore_dataset.R
# how this cutoff was chosen
filter.id = paste0("mincelllineWithData", min.cellline)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
rt <- readRDS(file=reptimePath)

# Rows of bins with data from cell lines < min.cellline are populated with NAs
rt.mx <- data.matrix(mcols(rt))
nonNAPerBin <- apply(rt.mx, MARGIN=1, FUN=function(bin){
  sum(!is.na(bin))
})
rt.mx[nonNAPerBin < min.cellline, ] <- NA

# Separate names of nontumor from tumor cell lines
cell.nmes <- dimnames(rt.mx)[[2]]
cell.nmes.split <- strsplit(cell.nmes, split="_")
disease <- unlist(lapply(cell.nmes.split, FUN=function(nmes){
  return( tail(nmes, n=1) )
}))

ntumor.ct <- cell.nmes[ disease == "NA" ]
tumor.ct <- setdiff(cell.nmes, ntumor.ct)

# Calculate mean, median values for all and seperately for nontumor and tumor cell lines
vals.mx <- matrix(data=NA, nrow=nrow(rt.mx), ncol=6,
                  dimnames=list(NULL,
                                c("mean.all", "median.all",  "mean.tumor", "median.tumor", 
                                  "mean.nontumor", "median.nontumor")))
  
vals.mx[,"mean.all"] <- rowMeans(rt.mx, na.rm=T)
vals.mx[,"median.all"] <- apply(X=rt.mx, MARGIN=1, FUN=median, na.rm=T)

vals.mx[,"mean.tumor"] <- rowMeans(rt.mx[,tumor.ct], na.rm=T)
vals.mx[,"median.tumor"] <- apply(X=rt.mx[,tumor.ct], MARGIN=1, FUN=median, na.rm=T)

vals.mx[,"mean.nontumor"] <- rowMeans(rt.mx[,ntumor.ct], na.rm=T)
vals.mx[,"median.nontumor"] <- apply(X=rt.mx[,ntumor.ct], MARGIN=1, FUN=median, na.rm=T)

# rowMeans(rt.mx, na.rm=T) gives NaN when taking mean with all NAs (for median, it returns NA)
vals.mx[is.nan(vals.mx)] <- NA

# Finalise format of rt data for downstream plotting
rt.transformed.df <- cbind.data.frame(chroms=seqnames(rt), starts=start(rt), vals.mx)

saveRDS(rt.transformed.df, file=paste0(out.dir, "/RT_data_hg19_", filter.id, ".RData"))

# rm(list=ls()); gc()