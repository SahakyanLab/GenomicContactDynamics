################################################################################
# Histogram of the distance between Cp=21 contacts (per chr)
# Mac, R/3.6.1
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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_Cp21_gapdist")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
#chr.v = paste0("chr", c(1:22, "X"), sep="")
# In order of decreasing number of Cp=21 contacts
#chr.v = c("chr2","chr4","chr5","chr3","chr6","chr8","chr13","chr1","chr9","chr7",
#          "chr11","chr12","chr10","chr14","chr18","chr20","chr15","chr16","chr21",
#          "chr17","chr19","chrX","chr22")
chr = "chr1"
HiC.res = 4e4L
prob = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
pdf(file=paste0(out.dir, "/", gcb, "_gapdist_Cp21.pdf"), width=30, height=15)
par(mfrow=c(4,6))

num.ij.lst <- list()

for(chr in chr.v){
  
  # Load PERSIST.MX
  load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  gap <- PERSIST.MX$hits[PERSIST.MX$ntis==21, c("i", "j")]
  rm(PERSIST.MX); gc()
  
  gap <- (gap$j-gap$i)*HiC.res/10^6
  
  num.ij <- length(gap)
  
  #write(paste(chr, num.ij, sep="\t"), append=T,
  #      file=paste0(out.dir, "/", gcb, "_gapdist_Cp21"))
  
  ylb <- ifelse(prob==TRUE, "Probability", "Frequency")
  # Main
  hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
       main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
       xlab=expression(bold("Distance between contacting regions, Mb")),
       ylab=bquote(bold( .(ylb) ))
  )
  #mtext(text=bquote(bold(.(num.ij)~"c"["p"]~"=21 contacts")),
  #      side=3, line=1e-4, at=4)
  
  print(paste0(chr, " done!"), quote=FALSE)
}

dev.off()

# rm(list=ls())