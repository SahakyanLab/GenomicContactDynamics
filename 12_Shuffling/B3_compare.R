################################################################################
# Compare discordance of original and shuffled HiC contacts
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
#whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    orig.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints/out_constraints/merged_final"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/8_ShuffleContactBins"
    orig.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints/out_constraints/merged_final"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
shuff.dir = paste0(wk.dir, "/out_constraints/merged_final")
out.dir = paste0(wk.dir, "/out_compare")
### OTHER SETTINGS #############################################################
chr.v = "chrALL" #paste("chr", c("ALL", 1:22, "X"), sep="")
gcb = "min2Mb"
kmer.len = 7L
type = "align" # kmer | align
affix = "_ijShuffled"
plotOnly = FALSE
mannwhit = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
### FUNCTION ###################################################################
HicHybridCompare <- function(
  orig.dir = "/dir",
  shuff.dir = "/dir",
  out.dir = "/dir",
  chr = "chr1",
  gcb = "min2Mb",
  kmer.len = 7L,
  # Identifier of the shuffled set
  affix = "_ijShuffled",
  plotOnly = FALSE,
  mannwhit = TRUE
){
  if(plotOnly==FALSE){
    set.id <- c("", affix)
    set.id.len <- length(set.id)
    # Make df with the following columns: Gfree-sdDifference-NegSumAbsDiff-cp-set
    # Set 1 = orig; Set 2 = shuffled
    HYBCOMB.DF <- list()
    for(i in 1:set.id.len){
      set <- set.id[i]
      dir <- ifelse(i==1, orig.dir, shuff.dir)
      # CII.MX
      load(file=paste0(dir, "/", chr, "_", type, "_", gcb, set, ".RData"))
      clnme <- colnames(CII.MX)
      # Remove contact without Cp
      CII.MX <- CII.MX[!is.na(CII.MX[,"Cp"]),]
      HYBCOMB.DF[[i]] <- data.frame(CII.MX[,!clnme%in%c("i", "j")], set=i)
      colnames(HYBCOMB.DF[[i]]) <- c(clnme[!clnme%in%c("i", "j")], "set")
      print(set, quote=FALSE)
      rm(CII.MX, set, dir); gc()
    }
    HYBCOMB.DF <- do.call("rbind", HYBCOMB.DF)
    save(HYBCOMB.DF, file=paste0(out.dir, "/", chr, "_HybComb_", type, "_", gcb, ".RData"))
  } else {
    load(paste0(out.dir, "/", chr, "_HybComb_", type, "_", gcb, ".RData"))
  }
  #---------------------------------------
  # Change C|| to CII to make it valid for the Mann-Whitney test
  clnme <- colnames(HYBCOMB.DF); clnme[clnme=="C||"] <- "CII"
  colnames(HYBCOMB.DF) <- clnme
  feat.v <- clnme[ !clnme%in%c("Cp", "set") ]
  if(type=="align"){
    ylab <- list(CII=bquote(bold( "c"["||"]~"align" )) 
                 )
  } else if(type=="kmer"){
    ylab <- list(CII=bquote(bold( "c"["||"]~"kmer" )),
                 Gfree=bquote(bold( "Gfree, kcal/mol" )),
                 sdDifference=bquote(bold( "s ("~"c"["||"]~")" ))
    )
  } else {
    stop("Invalid type.")
  }
  #---------------------------------------
  id <- paste0(chr, "_" , gcb, "_", type)
  cp.v <- sort( as.numeric(unique(HYBCOMB.DF$Cp)), decreasing=FALSE )
  cp.v.len <- length(cp.v)
  # Specify order to make sure that orig will always be on the left and shuff on the right
  HYBCOMB.DF$set <- factor(as.character(HYBCOMB.DF$set), levels=c("1", "2"))
  # Turn into a factor, whose levels are ordered from 1 to 21. This order determines the
  # sequence along the x-axis.
  HYBCOMB.DF$Cp <- factor(as.character(HYBCOMB.DF$Cp), levels=as.character(1:21))
  
  if(mannwhit){
    # Initialize output matrix with the p-values
    pval.mx <- matrix(data=NA, nrow=21, ncol=length(feat.v), dimnames=list(cp.v, feat.v))
  }
  #---------------------------------------
  if( any(is.na(HYBCOMB.DF$Cp)) ){
    stop("Checkpoint.")
  }
  # Count non-NA contacts orig and shuff to check the difference
  ij.orig.TF <- HYBCOMB.DF$set==1
  ij.shuff.TF <- !ij.orig.TF
  feat.nonNA.TF <- !is.na(HYBCOMB.DF$CII)
  #-------------------
  tot.ij.orig <- sum(ij.orig.TF)
  perc.ij.shuff <- sum(ij.shuff.TF)/tot.ij.orig*100
  perc.nonNA.ij.orig <- sum(ij.orig.TF & feat.nonNA.TF)/tot.ij.orig*100
  perc.nonNA.ij.shuff <- sum(ij.shuff.TF & feat.nonNA.TF)/tot.ij.orig*100
  #-------------------
  perCp.orig <- table(HYBCOMB.DF$Cp[feat.nonNA.TF & ij.orig.TF]); rm(ij.orig.TF)
  #percPerCp.orig <- perCp.orig/tot.ij.orig*100
  percPerCp.shuff <- (table(HYBCOMB.DF$Cp[feat.nonNA.TF & ij.shuff.TF])/perCp.orig)*100
  rm(ij.shuff.TF, feat.nonNA.TF)
  mx <- cbind(orig=c(tot=tot.ij.orig, nonNA=perc.nonNA.ij.orig, perCp.orig),
              shuff=c(tot=perc.ij.shuff, nonNA=perc.nonNA.ij.shuff, percPerCp.shuff))
  write.table(mx, file=paste0(out.dir, "/", id, "_counts"), col.names=TRUE, row.names=TRUE,
              quote=FALSE, sep="\t")
  rm(perCp.orig, percPerCp.shuff); gc()
  #---------------------------------------
  for(feat in feat.v){
    jpeg(filename=paste0(out.dir, "/", id, "_", feat, "_combBP.jpeg"),
         units="in", width=12, height=10, res=500)
    # default = par(mar=c(5,4,4,2)+0.1);  c(bottom, left, top, right)
    # default = par(mgp=c(3, 1, 0)); c(axis title, axis labels, axis line)
    par(mar=c(5.5, 6.5, 6.5, 2)+0.1, mgp=c(3, 1.5, 0)) 
    myplot <- boxplot(as.formula(paste0(feat,"~set*HYBCOMB.DF$Cp")), outline=FALSE,
                      data=HYBCOMB.DF, boxwex=0.6, xlab="", ylab="", main="", 
                      cex.axis=2.5, col=c("#FDC776" , "gray91"), xaxt="n")
    # Labelling x-axis
    # "1" "1" "2" "2" "3" "3"
    xlabs <- sapply( strsplit(x=myplot$names , split='\\.') , 
                     FUN=function(x) x[[2]] )
    # "1" "2" "3" "4"
    xlabs <- xlabs[seq(1, length(xlabs), 2)]
    axis(side=1, at=seq(1.5, cp.v.len*2 , by=2), 
         labels=1:21, tick=TRUE, cex.axis=2.5)
    
    # X-axis title
    mtext(side=1, text=expression( bold("c"["p"]) ), line=4.5, cex=3)
    # Y-axis title
    mtext(side=2, text=bquote(bold( .(ylab[[feat]]) )), line=3.5, cex=3)
    # Plot title
    mtext(side=3, text=paste0(id, "_outline=NA"), line=1.5, cex=2)
    
    # Line separating two boxplots per Cp
    for(i in seq(0.5 , cp.v.len*2+1, 2)){ abline(v=i, lty=1, col="gray50")}
    legend("topright", legend=c( expression(bold("orig")), expression(bold("shuff")) ), 
           col=c("#FDC776" , "gray91"), pch=15, bty="o", bg="white", pt.cex=3, 
           cex=2, horiz=F, inset=c(0,-0.17), xpd=TRUE)
    dev.off()
    #---------------------------------------
    if(mannwhit){
      # Wilcoxon-Mann-Whitney test (non-parametric, unpaired/independent groups) 
      # to compare orig and shuff values per Cp
      # The test assumes that the shape of the two distributions are similar so
      # check the boxplots
      # If both x and y are given and paired is FALSE, a Wilcoxon rank sum test
      # (equivalent to the Mann-Whitney test) is carried out.
    
      for(i in 1:cp.v.len){
        cp <- cp.v[i]
        mw <- wilcox.test(as.formula(paste0(feat, "~set")), alternative="two.sided",
                          data=HYBCOMB.DF[HYBCOMB.DF$Cp==cp & !is.na(HYBCOMB.DF[[feat]]), 
                                          c(feat, "set")],
                          paired=FALSE)
        # tt <- t.test("set", feat)
        pval.mx[as.character(i),colnames(pval.mx)==feat] <- mw$p.value
        rm(mw, cp)
      }
    }
    #---------------------------------------
    print(paste(chr, ":", feat, " done!"))
  } # feat.v for loop
  
  if(mannwhit){
    write.table(x=pval.mx, file=paste0(out.dir, "/", id, "_OlessS_mwtest"),
                col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
  }
}
################################################################################
HicHybridCompare <- cmpfun(HicHybridCompare, options=list(suppressUndefined=TRUE))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chr.v){
  HicHybridCompare(
    orig.dir=orig.dir,
    shuff.dir=shuff.dir,
    out.dir=out.dir,
    chr=chr,
    gcb=gcb,
    kmer.len=kmer.len, 
    affix=affix, 
    plotOnly=plotOnly,
    mannwhit=mannwhit
  )
}
# rm(list=ls())
  


