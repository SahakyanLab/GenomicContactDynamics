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
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
permtsumPath = paste0(wk.dir, "/out_summary/nperm10000_seed429_CptopCP3_permtsum.RData")
out.dir = paste0(wk.dir,"/out_table")
### OTHER SETTINGS #############################################################
out.name = "nperm10000_seed429_CptopCP3"
# Feature priority list name
pr.name = "nperm10000_seed429_CptopCP3"
pr.eval.v = c("numOlapA", "comOlap")
pval.cutoff = 0.05
regenerateCSV = TRUE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
pr.name <- paste0(pr.name, "_", paste(pr.eval.v, collapse="_"))
if(regenerateCSV){
  load(file=permtsumPath)
  # Initialise big table
  df <- as.data.frame(matrix(data=NA, nrow=nrow(PERMTSUM$pval), ncol=ncol(PERMTSUM$pval)*3))
  col.len <- ncol(df)
  df[,seq(from=1, to=col.len, by=3)] <- PERMTSUM$obs
  df[,seq(from=2, to=col.len, by=3)] <- PERMTSUM$alt
  df[,seq(from=3, to=col.len, by=3)] <- PERMTSUM$pval
  colnames(df) <- paste( rep(colnames(PERMTSUM$pval), each=3), 
                         rep(c("obs", "alt", "pval"), times=3), sep="_" )
  df <- cbind.data.frame(foi=PERMTSUM$foifile, df, stringsAsFactors=FALSE)
  rm(PERMTSUM); gc()
  write.csv(df, file=paste0(out.dir, "/", out.name, "_permtsum.csv"),
            row.names=FALSE)
} else {
  df <- read.csv(file=paste0(out.dir, "/", out.name, "_permtsum.csv"), header=TRUE,
                 check.names=FALSE, stringsAsFactors=FALSE)
}
# Identify priority features (significantly enriched based on total overlap and
# %number of overlapping cp regions)
alt.TF <- pval.TF <- rep(TRUE, times=nrow(df))
for(pr in pr.eval.v){
  alt.TF <- alt.TF & df[,paste0(pr, "_alt")]=="greater" 
  pval.TF <- pval.TF & df[,paste0(pr, "_pval")] < pval.cutoff
  print(pr, quote=FALSE)
}
#pr.foi.v <- df$foi[alt.TF & pval.TF]
writeLines(text=df$foi[alt.TF & pval.TF], 
           con=paste0(out.dir, "/foifile_priority_", pr.name))

#temp <- unlist(
#  lapply(X=strsplit(x=pr.foi.v, split="_foi_|\\_"), 
#         FUN=function(f){f[3]})
#)

#pr.foi.v <- pr.foi.v[duplicated(temp)]
#writeLines(text=pr.foi.v, con=paste0(out.dir, "/foifile_priority_", pr.name, "_acrossTissues"))

# rm(list=ls()); gc()



