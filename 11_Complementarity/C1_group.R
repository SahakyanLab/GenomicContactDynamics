################################################################################
# Categorise complementarity values into 3 categories (high, middle, low)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# 
data.dir = out.dir = paste0(wk.dir, "/out_constraints/merged_final")
out.dir = paste0(wk.dir, "/out_group")
### OTHER SETTINGS #############################################################
chr.v = paste("chr", c(1:22, "X"), sep="")
gcb = "min2Mb"
type = "align"
cutoff = 5 # percentage
writeCsv = TRUE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chr.v){
  # Load CII.MX
  load(file=paste0(data.dir, "/", chr, "_", type, "_", gcb, ".RData"))
  val <- CII.MX[!is.na(CII.MX[,"C||"]),"C||"]
  
  pdf(file=paste0(out.dir, "/", chr, "_", type, "_", gcb, "_grouped_cutoff", cutoff, ".pdf"),
      width=10, height=10)
  par(mfrow=c(2,2))
  
  len <- length(val)
  tot.ij <- length(CII.MX[,1])
  nonNA.ij <- format(len/tot.ij*100, scientific=TRUE, digits=4)
  tot.ij <- format(as.numeric(length(CII.MX[,1])), scientific=TRUE, digits=4)
  
  hist(val, col="#55bde6", xlab=expression("C"["||"]), breaks=30,
       main=paste0(chr, "_", gcb, "_", tot.ij, "ij_", nonNA.ij, "nonNACII"))
  #---------------------------------------
  # Categorize continuous variable, C||, into 1 (top 5%), -1 (bottom 5%),
  # 0 (rest) - based on area under the curve (frequency)
  # Note that the 5% cut-off value was determined using the total number of contacts
  # with C|| value (meaning excluding those with missing C|| values)
  # Ideal 5% cut-off based on the number of datapoints
  val.srt <- sort(val, decreasing=TRUE)
  rm(val, tot.ij, nonNA.ij); gc()
  
  cutval <- ceiling(len*cutoff/100) # = 782741
  
  # Upper half minimum based on ideal cut-off
  up.b <- val.srt[cutval] # = -0.9028
  # unique(val.srt): -0.9026 -0.9026 -0.9027 -0.9028 -0.9029 -0.9030
  
  # Lower half minimum based on ideal cut-off
  low.b <- val.srt[len-cutval] # = -1.7029
  # unique(val.srt): -1.7028 -1.7028 -1.7029 -1.7029 -1.7030 -1.7030
  
  # Because many contacts have the same C|| value, we just use the minimums based
  # on the ideal cut-off for the grouping instead of using the actual cut-off = 782741
  group <- CII.MX[,4]
  # N = 781810
  group[CII.MX[,4] > (up.b)] <- 1
  # N = 14090461
  group[(CII.MX[,4] <= (up.b)) & (CII.MX[,4] >= (low.b))] <- 0
  # N = 782539
  group[CII.MX[,4] < low.b] <- -1
  #---------------------------------------
  CII.MX <- cbind(CII.MX, group)
  rm(cutval, up.b, low.b, group, val.srt); gc()
  
  # Check, should be TRUE
  identical( which(is.na(CII.MX[,"C||"])), which(is.na(CII.MX[,"group"])) )
  
  # Histograms
  hist( CII.MX[CII.MX[,"group"]==1,"C||"], main=paste0("group 1_", 
                                                       sum(CII.MX[,"group"]==1, na.rm=TRUE),
                                                       "_cutoff=", cutoff, "%"), 
        xlab=expression("C"["||"]), col="#55bde6")
  hist( CII.MX[CII.MX[,"group"]==0,"C||"], main=paste0("group 0_",
                                                       sum(CII.MX[,"group"]==0, na.rm=TRUE),
                                                       "_cutoff=", cutoff, "%"),
        xlab=expression("C"["||"]), col="#55bde6")
  hist( CII.MX[CII.MX[,"group"]==-1,"C||"], main=paste0("group -1_",
                                                        sum(CII.MX[,"group"]==-1, na.rm=TRUE),
                                                        "_cutoff=", cutoff, "%"),
        xlab=expression("C"["||"]), col="#55bde6")
  dev.off()
  
  save(CII.MX, file=paste0(out.dir, "/", chr, "_", type, "_", gcb, "_grouped_cutoff", cutoff, ".RData"))
  
  if(writeCsv){
    write.csv(CII.MX, file=paste0(out.dir, "/", chr, "_", type, "_", gcb, "_grouped_cutoff", cutoff, ".csv"),
              quote=FALSE, row.names=FALSE)
  }
  
  rm(CII.MX); gc()
  
  print(paste0(chr, " done!"), quote=FALSE)
}

#x <- read.csv(file=paste0(out.dir, "/", chr, "_", type, "_", gcb, "_grouped.csv"))
#identical( which(is.na(x[,"C.."])), which(is.na(x[,"group"])) )
#NAval <- x[is.na(x[,"group"]), 4]
## Check, should be true
#all(is.na(NAval))

# rm(list=ls())
