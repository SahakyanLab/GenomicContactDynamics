################################################################################
# Make a csv file with the values from PERMTSUM
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    binmx.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Circos/out_bindata"
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
permtsumPath = paste0(wk.dir, "/out_summary_runA_cp21_b1b2b3/nperm10000_cp21_seed342_runA_cp21_b1b2b3_permtsum.RData")
out.dir = paste0(wk.dir,"/out_table")
### OTHER SETTINGS #############################################################
out.name = "nperm10000_cp21_seed342_runA_cp21_b1b2b3"
pr.name = "runA_numOlapBANDcomOlap"
priority = c("% number of B regions overlapping", "total length of intersection")
pval.cutoff = 0.05
regenerateCSV = TRUE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(regenerateCSV){
  load(file=permtsumPath)
  # Initialise big table
  df <- as.data.frame(matrix(data=NA, nrow=nrow(PERMTSUM$pval), ncol=ncol(PERMTSUM$pval)*3))
  
  df[,seq(from=1, to=ncol(df), by=3)] <- PERMTSUM$obs
  df[,seq(from=2, to=ncol(df), by=3)] <- PERMTSUM$alt
  df[,seq(from=3, to=ncol(df), by=3)] <- PERMTSUM$pval
  colnames(df) <- paste( rep(colnames(PERMTSUM$pval), each=3), 
                         rep(c("obs", "alt", "pval"), times=3), sep="_" )
  df <- cbind.data.frame(foi=PERMTSUM$foifile, df, stringsAsFactors=FALSE)
  rm(PERMTSUM)
  write.csv(df, file=paste0(out.dir, "/", out.name, "_permtsum.csv"),
            row.names=FALSE)
} else {
  df <- read.csv(file=paste0(out.dir, "/", out.name, "_permtsum.csv"), header=TRUE,
                 check.names=FALSE, stringsAsFactors=FALSE)
}

# Identify priority features (significantly enriched based on total overlap and
# %number of overlapping cp regions)
alt.TF <- pval.TF <- rep(TRUE, times=nrow(df))
for(pr in priority){
  alt.TF <- alt.TF & df[,paste0(pr, "_alt")]=="greater" 
  pval.TF <- pval.TF & df[,paste0(pr, "_pval")] < pval.cutoff
  print(pr, quote=FALSE)
}

writeLines(text=df$foi[alt.TF & pval.TF], con=paste0(out.dir, "/foifile_priority_", pr.name))

# rm(list=ls())



