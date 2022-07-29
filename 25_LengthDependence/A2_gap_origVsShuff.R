################################################################################
# Compare contact gap distributions of orig (real) and shuffled contacts
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
    home.dir = "/project/sahakyanlab/ltamon" 
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/25_LengthDependence")
persistReal.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
persistShuf.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/12_Shuffling/out_features")
out.dir = paste0(wk.dir, "/out_gap_origVsShuff")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
bin.len = 40000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chr.v){
  
  load(paste0(persistReal.dir, "/", chr, "_Persist_", gcb, ".RData"))
  df <- cbind.data.frame(gap=PERSIST.MX$hits$j - PERSIST.MX$hits$i, 
                         Cp=PERSIST.MX$ntis, set="orig", stringsAsFactors=F)
  rm(PERSIST.MX)
  
  #
  load(paste0(persistShuf.dir, "/", chr, "_Persist_", gcb, "_ijShuffled.RData"))
  df <- rbind.data.frame(df,
                         cbind.data.frame(gap=PERSIST.MX$hits[,"j"] - PERSIST.MX$hits[,"i"] - 1, 
                                          Cp=PERSIST.MX$ntis, set="shuff", stringsAsFactors=F))
  rm(PERSIST.MX); gc()
  
  df$set <- factor(x=as.character(df$set), levels=c("orig", "shuff"))
  
  #
  save(df, file=paste0(out.dir, "/", chr, "_", gcb, "_", bin.len, "_gap_real_shuff_tmp.RData"))
  rm(df); gc()
  
  print(paste0(chr, " done!"), quote=F)
  
}

df <- sapply(X=chr.v, simplify=F, FUN=function(chr){
  df.file <- paste0(out.dir, "/", chr, "_", gcb, "_", bin.len, "_gap_real_shuff_tmp.RData")
  load(df.file)
  file.remove(df.file)
  return(df)
})
df <- do.call("rbind.data.frame", df)

out.id <- paste0("chrALL_", gcb, "_binlength", bin.len, "bp")

jpeg(filename=paste0(out.dir, "/", out.id, "_gap_real_shuff_combBP.jpeg"),
     units="in", width=12, height=10, res=500)
boxplot(formula=gap~set*Cp, outline=F, data=df, boxwex=0.6, xlab="Cp", cex.axis=0.2,
        ylab="Contact gap j-i-1, bins", main=paste0(out.id, "_Nij=", length(df$Cp)), 
        col=c("#FDC776" , "gray91"))
legend("topright", legend=c( expression(bold("orig")), expression(bold("shuff")) ), 
       col=c("#FDC776" , "gray91"), pch=15, bty="n", bg="white", pt.cex=3, 
       cex=2, horiz=F, inset=c(0,0), xpd=T)
dev.off()

# rm(list=ls()); gc()