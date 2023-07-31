################################################################################
# Prediction vs. truth
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/26_PredictingCp")
persist.dir  = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
mltbl.dir = paste0(wk.dir, "/z_ignore_git/out_makeTable")
pred.dir = paste0(wk.dir, "/out_autogluon_run") # same order as MLtbl
out.dir  = paste0(wk.dir, "/out_assessPrediction") # Path to output directory 
### OTHER SETTINGS #############################################################
model_id = "agModel_min2Mb_chr1_grpcomplgrpkmer3grpkmer1grpGC_best_quality_18000sec_regression_subsample_seed234_5percExceptgrthn18_1237618datapoints"
gcb = "min2Mb"
chrs_test = "chr18"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Get pre

preds <- fread(file=paste0(pred.dir, "/", model_id, "_chrtest_", chrs_test, "_MLtbl_pred.csv"),
               stringsAsFactors=F, data.table=F, header=T)

truths <- fread(file=paste0(mltbl.dir, "/", gcb, "_", chrs_test, "_MLtbl.csv"),
                stringsAsFactors=F, data.table=F, header=T)$LABEL

if( !identical(0:(length(truths)-1), preds$V1) ){
  rm(preds); stop("preds and truths not same order.")
} 

preds <- preds$LABEL
truths <- factor(as.character(truths), levels=as.character(1:21))

boxplot(formula=preds~truths, data=cbind.data.frame(preds=preds, truths=truths))

# Perform different correlation tests to serve as performance metric

cor.test(x=as.numeric(as.character(truths)), y=preds, alternative="two.sided", method="spearman")

# rm(list=ls()); gc()