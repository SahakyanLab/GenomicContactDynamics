################################################################################
# Boxplot of the contact gap lengths for all long-range contacts (Cp=0) and
# per Cp. Plot combines contacts from all chr. The boxplot whiskers is set to 
# include min and max value so there should be no outliers. Boxplot statistics are 
# saved in the csv file. Plot combines contacts from all chr. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/3_TADboundary"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster") {
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/3_TADboundary"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# Bed files directory
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_gapdist")
### OTHER SETTINGS #############################################################
gcb = "min2Mb" #"min05Mb"
chr.v = "chr21" #paste0("chr", c(1:22, "X"))
Cp.v = 1:21
bin.len = 4e4
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
### FUNCTION ###################################################################
makebp <- function(out.name=paste0(out.dir, "/", out.name, "_AllLR"),
                   GAP=GAP, xlab=levels(GAP$Cp)
){
  xlab.len <- length(xlab)
  max.y <- ceiling(max(GAP$gap))
  
  png(file=paste0(out.name, ".png"), width=10, height=10, units="in", res=1200)
  # range=0 causes the whiskers to extend to the data extremes (and no
  # outliers be returned)
  bp <- boxplot(gap~Cp, outline=TRUE, data=GAP,  xlab=expression(bold("c"["p"])), 
                ylab=expression(bold("Mb")), width=rep(0.4, times=xlab.len), 
                cex.axis=1, col="#FDC776", xaxt="n", range=0, ylim=c(0,max.y))
  rm(GAP);
  axis(side=1, at=1:xlab.len, labels=xlab, cex.axis=1)
  mtext(side=3, line=2, cex=0.7, 
        text=paste0(out.name, 
                    "\n whiskers stretch to minmax (range=0); coef not set but should be 1.5(overrode by range=0?);
                    Cp=0 includes all LR contacts; gapdist does not include contact regions"))
  dev.off()
  
  return(bp$stats)
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, "..."), quote=FALSE)

Cp.v <- sort(unique(Cp.v), decreasing=FALSE)
Cp.v.len <- length(Cp.v)
GAP <- vector(mode="list", length=length(Cp.v))
names(GAP) <- Cp.v

for(chr in chr.v){
  
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  gap.v <- (PERSIST.MX$hits$j-PERSIST.MX$hits$i-1L)
  cp.v <- PERSIST.MX$ntis
  rm(PERSIST.MX)
  gc()
  
  for(Cp in Cp.v){
    GAP[[as.character(Cp)]] <- c(GAP[[as.character(Cp)]], gap.v[cp.v==Cp])
  }
  
  rm(gap.v, cp.v); gc()
  print(paste0(chr, " done!"), quote=FALSE)
  
}

GAP <- stack(GAP)
colnames(GAP) <- c("gap", "Cp")
GAP$gap <- GAP$gap*(bin.len/1e6)
out.name <- paste0(gcb, "_LR_gapdist_bp")

#-------------------Boxplot across Cp

GAP$Cp <- factor(as.character(GAP$Cp), levels=c("0",as.character(Cp.v)))
bp.stat <- makebp(out.name=paste0(out.dir, "/", out.name, "_acrossCp"), GAP=GAP, 
                  xlab=levels(GAP$Cp))
bp.stat <- bp.stat[,-1] # Remove stat for factor "0", which are all NAs
colnames(bp.stat) <- as.character(Cp.v)
rownames(bp.stat) <- c("Min", "25th", "Median", "75th", "Max")

#-------------------Boxplot for all LR contacts

GAP$Cp[GAP$Cp%in%Cp.v] <- "0"
temp <- makebp(out.name=paste0(out.dir, "/", out.name, "_All"), GAP=GAP, 
               xlab=levels(GAP$Cp))
bp.stat <- cbind(AllLR=temp[,1], bp.stat)
rm(temp, GAP)
write.csv(x=bp.stat, file=paste0(out.dir, "/", out.name, "_stat_Mbunit.csv"), 
          row.names=TRUE, quote=FALSE)

# rm(list=ls()); gc()
