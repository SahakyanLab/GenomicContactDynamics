################################################################################
# Combine shuffled contact pairs per Cp into one datafile with format same as
# PERSIST.MX. NOTE: Numbers below QEO plot is the number of contacts per Cp in 
# the original set
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    persist.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else if(whorunsit == "LiezelLinuxDesk"){
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    persist.dir = "/home/ltamon/Database/HiC_features_GSE87112_RAWpc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(wk.dir, "/out_shuffled")
out.dir = paste0(wk.dir, "/out_features")
### OTHER SETTINGS #############################################################
# chr 1 min 2Mb cp12:21 - 3.723G - 46s
chr.v = paste("chr", c(1:22, "X"), sep="")
gcb = "min2Mb"
cp.v = 1:21 #1:21
affix = "_ijShuffled"
plotOnly = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
### FUNCTION ###################################################################
makeShuffledPersistMx <- function(
  persist.dir=persist.dir,
  shuffled.dir=data.dir, 
  out.dir=out.dir,
  gcb = "min2Mb",
  chr = "chr1", 
  affix = "_ijShuffled",
  cp.v = 1:21,
  plotOnly = FALSE
){
  
  print(paste0(chr, "..."), quote=FALSE)
  
  id <- paste0(gcb, "_", chr)

  if(plotOnly==FALSE){
    
    # Required minimum gap between contacting bins
    # In terms of bins
    mingap <- ifelse(gcb=="min2Mb", 50, 12.5)
    
    # Load PERSIST.MX 
    load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    ij.orig <- data.frame(PERSIST.MX$hits[,c("i", "j")], cp=PERSIST.MX$ntis, 
                          row.names=NULL, stringsAsFactors=FALSE)
    rm(PERSIST.MX); gc()
    
    QEO <- list()
    PERSIST.MX <- list()
    for(cp in cp.v){
      
      print(paste0("cp=", cp))
      
      #---------------------------------------------------------------------------
      # Picking the best replica
      #---------------------------------------------------------------------------
      reps.nme <- list.files(path=data.dir, pattern=paste0(id, "_cp", cp, "_"), 
                             ignore.case=TRUE, include.dirs=FALSE)
      reps.nme <- reps.nme[ grepl(x=reps.nme, pattern="O\\.Rdata") ]
      # O=-Q=E; O is the error rate
      O.stored.v <- sapply(X=reps.nme, FUN=function(rep){
        # O.stored from Optimus (% of contacts not satisfying the requirements)
        load(paste0(data.dir, "/", rep))
        return(O.stored)
      })
      bestrep.ind <- which(O.stored.v==min(O.stored.v))[1] 
      bestrep.O <- O.stored.v[bestrep.ind]
      names(bestrep.O) <- NULL
      bestrep.nme <- gsub(x=reps.nme[bestrep.ind], pattern="O", replacement="K",
                          fixed=TRUE)
      print(bestrep.nme)
      rm(O.stored.v, bestrep.ind)
      
      # Load original ij pairing
      K.original <- ij.orig[ij.orig$cp==cp, c("i", "j")]
      len <- nrow(K.original)
      K.original <- list(hits=cbind(i=K.original$i, j=K.original$j))
      # Best shuffled ij pairs (K.stored)
      load(paste0(data.dir, "/", bestrep.nme))
      # Check
      if( !identical(K.original$hits[,"i"], K.stored$hits[,"i"]) ){
        stop("i's of original and shuffled not identical.")
      } 
      # Check
      if(len!=K.stored$len){
        stop("Lengths of original and shuffled not equal.")
      } 
      # Manual calculation of O for comparison with O calculated by Optimus
      mg.log <- (K.stored$hits[,"j"]-K.stored$hits[,"i"]) <= mingap
      log.init <- mg.log | ( K.stored$hits[,"j"]==K.original$hits[,"j"] )
      O.opti <- (sum(log.init)/len)*100
      
      if(O.opti!=bestrep.O){
        stop("O.opti not equal to O.stored!")
      } 
      
      # Update O.opti accounting for duplicates
      log.init <- log.init | duplicated(K.stored$hits)
      O.opti <- ( sum(log.init)/len )*100
      
      # Attemp to reduce errors in best shuffled data from Optimus
      K.temp <- K.stored
      # Reverse pairs with i > j
      K.temp$hits[mg.log,] <- cbind(K.stored$hits[mg.log, "j"], K.stored$hits[mg.log, "i"])
      rm(mg.log); gc()
      # Calculate new O
      K.orig.str <- paste(K.original$hits[,"i"], K.original$hits[,"j"], sep="_")
      K.temp.str <- paste(K.temp$hits[,"i"], K.temp$hits[,"j"], sep="_")
      ind <- match(x=K.temp.str, table=K.orig.str)
      # duplicated contacts | present in original set | below mingap
      log.temp <- duplicated(K.temp.str) | !is.na(ind) | (K.temp$hits[,"j"]-K.temp$hits[,"i"]) <= mingap 
      O.temp <- (sum(log.temp)/len)*100
      
      if(O.temp < O.opti){
        print("K.temp used")
        O.fin <- O.temp
        log.fin <- log.temp
        K.stored <- K.temp; rm(K.temp, O.temp, log.temp); gc()
      } else {
        print("K.stored used")
        log.fin <- log.init 
        O.fin <- ( sum(log.fin)/len )*100
        rm(K.temp, O.temp, log.init); gc()
      }
      print(paste0("%Error=", O.fin))
      rm(K.original); gc()
      QEO[[as.character(cp)]] <- c(cp, len, O.opti, O.fin)
      # Remove contacts not satisfying requirements
      PERSIST.MX[[as.character(cp)]] <- cbind(K.stored$hits[!log.fin,], ntis=rep(cp))
      rm(K.stored, log.fin); gc()
      
    } # cp.v for loop end
    
    PERSIST.MX <- do.call("rbind", PERSIST.MX)
    PERSIST.MX <- list( hits=PERSIST.MX[,c("i", "j")], ntis=PERSIST.MX[,"ntis"] )
    
    save(PERSIST.MX, file=paste0(out.dir, "/", chr, "_Persist_", gcb, affix, ".RData"))
    rm(PERSIST.MX); gc()
    
    QEO <- do.call("rbind", QEO)
    colnames(QEO) <- c("cp", "N", "O.opti", "O.fin")
    
    save(QEO, file=paste0(out.dir, "/", chr, "_QEO_", gcb, affix, ".RData"))
    
  } else {
    load(file=paste0(out.dir, "/", chr, "_QEO_", gcb, affix, ".RData"))
  }
  
  # Plot O.opti and O.fin Vs Cp per chr
  #pdf(file=paste0(out.dir, "/", id, affix, "_QEOplot.pdf"), 
  #    width=10, height=10)
  #par(oma=c(0,0,0,0), mar=c(4.7, 5, 5, 2) + 0.1, mgp=c(3, 0.9, 0))
  plot(x=QEO[,"cp"], y=QEO[,"O.opti"], main="", cex.axis=2.2, xaxt="n",
       xlab="", ylab="", type="p", pch=19, cex=2, col="gray80",
       ylim=c( 0, ceiling( max(QEO[,-(1:2)]) ) )
       #ylim=c( floor( min(QEO[,-(1:2)]) ), ceiling( max(QEO[,-(1:2)]) ) )
  )
  points(x=QEO[,"cp"], y=QEO[,"O.fin"], pch=19, cex=2, col="black")
  axis(side=1, at=1:21, labels=1:21, cex.axis=2.2, mgp=c(3, 1.4, 0))
  # X axis title
  mtext(side=1, text=expression(bold( "c"["p"] )), line=3.4, cex=2)
  # Y axis title
  mtext(side=2, text=expression(bold( "% Error" )), line=3, cex=2)
  # Plot title
  mtext(side=3, text=paste0(id, affix, "_O.opti(with dupl)"), 
        line=2, cex=1.5)
  text(x=cp.v, y=0, labels=as.character(QEO[,"N"]), cex=0.5)
  #legend("topright", legend=c("Opti", "Fin"), inset=.02,
  #       col=c("gray80", "black"), pch=19, cex=1)
  legend("topright", legend=c("Opti", "Fin"), col=c("gray80", "black"),
         pch=19, bty="o", bg="white", pt.cex=2, cex=1.5, horiz=F, 
         # inset=c(0, -0.16) if 10x8
         inset=c(0, -0.12), xpd=TRUE)
  #dev.off()
  
  print(paste0(chr, " done!"), quote=FALSE)
  
}
################################################################################
makeShuffledPersistMx <- cmpfun(makeShuffledPersistMx, options=list(suppressUndefined=TRUE))
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
pdf(file=paste0(out.dir, "/", gcb, affix, "_QEOplot.pdf"), width=60, height=40)
par(mfrow=c(4,6))
for(chr in chr.v){
  makeShuffledPersistMx(
    persist.dir=persist.dir,
    shuffled.dir=data.dir, 
    out.dir=out.dir,
    gcb=gcb,
    chr=chr, 
    affix=affix,
    cp.v=cp.v,
    plotOnly=plotOnly
  )
}
dev.off()
# rm(list=ls())

