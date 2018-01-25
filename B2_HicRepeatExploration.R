## FUNCTION ####################################################################
HicRepeatExploration <- function(
  # Number of CPU cores to be used
  nCPU = 4, # optimal for the task on bucephalus (takes ~5GB per core for chr1)
  # Chromosome for which the data are to be pulled:
  chr = "chr1",
  # Feature database path and filename suffix:
  featureDBPath = "/home/alex/Desktop/CHROMSEQ/OUT",
        #"/Volumes/Data/Database/HiC_features_GSE87112_RAWpc",
  # Feature database filename suffix
  suffix = "min2Mb",
  # Wheteher to do the initial chr related expl. plots (contact vs. tissue freq)
  initialExplPlots = TRUE,
  # Whether to do the rest of the mobile DNA exploration
  mobDNAExplPlots = TRUE
){
################################################################################
library(hexbin)
library(RColorBrewer)
library(doMC)
library(foreach)
library(itertools)
registerDoMC(cores=nCPU)
################################################################################

id <- paste0(chr,"_",suffix)

###load(paste0(featureDBPath,"/",chr,"_BinKmer4_",suffix,".RData"))
###load(paste0(featureDBPath,"/",chr,"_BinKmer7_",suffix,".RData"))
load(paste0(featureDBPath,"/",chr,"_BinRep_",suffix,".RData"))
load(paste0(featureDBPath,"/",chr,"_Persist_",suffix,".RData"))

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

###########################
if(initialExplPlots==TRUE){

  # As of 25Jan18, hexbinplot does not work with pdf/jpeg/png devices while
  # both are specified inside a function. Hence, replotting would be needed
  # from PERSIST.MX in featureDBPath.

  df <- data.frame(ntis=PERSIST.MX$ntis, valsum=PERSIST.MX$valsum)
  pdf(file=paste0(id,"_hexplotA.pdf"), width=7, height=7)
  hexbinplot(valsum~ntis, data=df, colramp=rf,
             trans=log, inv=exp, aspect=1, mincnt=1,
             xlab="Num(non-0 tissues/cell lines)",
             ylab="Sum",
             main=id)
  dev.off()
  rm(df)

  df <- data.frame(ntis=PERSIST.MX$ntis, meanvalsum=PERSIST.MX$valsum/PERSIST.MX$ntis)
  pdf(file=paste0(id,"_hexplotB.pdf"), width=7, height=7)
  hexbinplot(meanvalsum~ntis, data=df, colramp=rf,
             trans=log, inv=exp, aspect=1, mincnt=1,
             xlab="Num(non-0 tissues/cell lines)",
             ylab="Sum / N(non-0)",
             main=id)
  dev.off()
  rm(df)

  df <- data.frame(ntis=PERSIST.MX$ntis, meanvalsum=PERSIST.MX$valsum/PERSIST.MX$ntis)
  pdf(file=paste0(id,"_boxplotB.pdf"), width=9, height=7)
  boxplot(formula=meanvalsum~ntis, data=df, col=rf(max(PERSIST.MX$ntis)),
          cex=0.1, ylim=c(1,14),
          outline=FALSE, xlab="Num(non-0 tissues/cell lines)",
          ylab="Sum / N(non-0)",
          main=id)
  dev.off()
  rm(df)

}
###########################

############################
if(mobDNAExplPlots == TRUE){

  element.names <- names(BINREP.MX[1,4:58])
  hits.len <- length(PERSIST.MX$hits[,1])
  ij.mx <- cbind(i=PERSIST.MX$hits[,"i"], j=PERSIST.MX$hits[,"j"])
  ntis.vec <- PERSIST.MX$ntis

  rm(PERSIST.MX)
  gc()

  ###---------------------------
  for(element in element.names){

    print(element, quote=FALSE)

    #### FOREACH EXECUTION #########
    element.count <- foreach(itr=isplitVector(1:hits.len, chunks=nCPU),
                             .combine="c", .inorder=TRUE,
                             .export=c("BINREP.MX","ij.mx","element"),
                       .noexport=ls()[!ls()%in%c("BINREP.MX","ij.mx","element")]
                            ) %dopar% {

       element.count.chunk <- sapply(itr,
                                     FUN=function(itr){
                                       min(
                                           BINREP.MX[which(BINREP.MX[,"bins"] %in%
                                                             ij.mx[itr,]),element]
                                          )
                                      },
                                      simplify=TRUE, USE.NAMES=FALSE)
      return(element.count.chunk)
    }
    ### END OF FOREACH EXECUTION ###

    df <- data.frame(ntis=ntis.vec, minELM=element.count)
    pdf(file=paste0(id,"_boxplot_",element,".pdf"), width=9, height=7)
      boxplot(formula=minELM~ntis, data=df,
              col=rf(max( ntis.vec )), cex=0.1, # ylim=c(0,30)
              outline=FALSE, xlab="Num(non-0 tissues/cell lines)",
              ylab="min repeat count in contact pairs",
              main=paste0(id,"_",element))
    dev.off()

    unique.ntis <- unique(ntis.vec)
    unique.ntis <- unique.ntis[order(unique.ntis)]
    num.contacts.each.ntis <- rep(NA,unique.ntis)
    num.non0cont.each.ntis <- rep(NA,unique.ntis)
    for(k in 1:length(unique.ntis)){
      num.contacts.each.ntis[k] <- sum( ntis.vec==unique.ntis[k] )

      num.non0cont.each.ntis[k] <- sum( (ntis.vec==unique.ntis[k]) &
                                        (element.count>=1) )
    }
    pdf(file=paste0(id,"_lineplot_",element,".pdf"), width=9, height=7)
      plot(x=unique.ntis, y=100*num.non0cont.each.ntis/num.contacts.each.ntis,
           col=rf(max(unique.ntis)), cex=1, pch=21, lwd=3, # ylim=c(0,30)
           xlab="Num(non-0 tissues/cell lines)",
           ylab="contact % with >=1 repeat pair",
           main=paste0(id,"_",element))
      lines(x=unique.ntis, y=100*num.non0cont.each.ntis/num.contacts.each.ntis)
    dev.off()

    # Dumped per element, then collected later, for memory efficiency.
    save(element.count, file=paste0(featureDBPath,"/",chr,
                                    "_MinElmSOURCE_",suffix,"_",element,".RData"))

    rm(df,
       element.count,
       unique.ntis,
       num.contacts.each.ntis,
       num.non0cont.each.ntis, k); gc()

  }
  ###---------------------------

  rm(ij.mx); gc()

  # Saving all data for min(repeat counts) in each interacting pairs of chr loci
  MINELM.MX <- matrix( NA, nrow=hits.len, ncol=length(element.names)+1 )
  dimnames(MINELM.MX)[[2]] <- c("ntis", element.names)
  MINELM.MX[,"ntis"] <- ntis.vec; rm(ntis.vec)
  for(element in element.names){
    SRC<-paste0(featureDBPath,"/",chr,"_MinElmSOURCE_",suffix,"_",element,".RData")
    load(SRC); file.remove(SRC)
    MINELM.MX[,element] <- element.count
    rm(element.count)
    gc()
  }
  save(MINELM.MX, paste0(featureDBPath,"/",chr,"_MinElm_",suffix,".RData"))

}
############################

print("HicRepeatExploration is DONE!", quote=FALSE)

}
################################################################################
suppressWarnings(suppressPackageStartupMessages(library(compiler)))
HicRepeatExploration <- cmpfun(HicRepeatExploration, options=list(suppressUndefined=TRUE))


# Execution of the function writen above...
################################################################################
HicRepeatExploration(
  # Number of CPU cores to be used
  nCPU = 4,
  # Chromosome for which the data are to be pulled:
  chr = "chr1",
  # Feature database path and filename suffix:
  featureDBPath = "/home/alex/Desktop/CHROMSEQ/OUT",
        #"/Volumes/Data/Database/HiC_features_GSE87112_RAWpc",
  # Feature database filename suffix
  suffix = "min2Mb",
  initialExplPlots = TRUE,
  mobDNAExplPlots = TRUE
)


















#library(ggplot2)
#p <- ggplot(df, aes(ntis,valsum)) +
#      stat_bin2d(binwidth=c(1,10)) +
#      scale_fill_gradientn(colours=r)



#> names(PERSIST.MX)
#[1] "hits"    "ntis"    "valsum"  "control"

#> PERSIST.MX$hits[1,]
#       i  j Co Hi Lu LV RV Ao PM Pa Sp Li SB AG Ov Bl MesC MSC NPC TLC ESC FC LC
#179619 2 54  0  0  0  0  0  0  0  0  0  0  0  0  0  0    0   1   0   0   0  1  0

#> PERSIST.MX$ntis
#[1]  2  7  2 10 14  7 12 ...

#> PERSIST.MX$valsum
#[1]   2  14   2  23  35  ...

#> PERSIST.MX$control[1,]
#i  j value
#172840 1 52     0



#> BINKMER4.MX[1,]
#bins startpos   endpos     AAAA     AAAC     AAAG     AAAT
#2    40001    80000       NA       NA       NA       NA   # NA means there is at least one N

#> BINREP.MX[1,]
#bins  startpos  endpos    ERV1    L1  ...
#2          40001          80000   0   ...

#> names(BINREP.MX[1,])
#[1] "bins"           "startpos"       "endpos"         "ERV1"           "L1"
#[6] "Alu"            "ERVL"           "hAT-Charlie"    "ERVL-MaLR"      "MIR"
#[11] "Satellite"      "ERVK"           "Low_complexity" "L2"             "TcMar-Tigger"
#[16] "Simple_repeat"  "hAT-Blackjack"  "Unknown"        "centr"          "RTE"
#[21] "TcMar-Mariner"  "Gypsy?"         "CR1"            "LTR"            "hAT-Tip100"
#[26] "hAT"            "Gypsy"          "rRNA"           "DNA"            "tRNA"
#[31] "TcMar?"         "srpRNA"         "TcMar-Tc2"      "ERV"            "snRNA"
#[36] "Helitron"       "hAT?"           "RTE-BovB"       "SINE"           "ERVL?"
#[41] "PiggyBac"       "LTR?"           "DNA?"           "Other"          "Merlin"
#[46] "MuDR"           "SINE?"          "Deu"            "RNA"            "scRNA"
#[51] "TcMar"          "Dong-R4"        "PiggyBac?"      "Unknown?"       "Penelope?"
#[56] "L1?"            "Helitron?"      "telo"

#
# [9] "chr10_BinKmer4_min05Mb.RData" "chr10_BinKmer4_min2Mb.RData"
# [11] "chr10_BinKmer7_min05Mb.RData" "chr10_BinKmer7_min2Mb.RData"
# [13] "chr10_BinRep_min05Mb.RData"   "chr10_BinRep_min2Mb.RData"
# [15] "chr10_Persist_min05Mb.RData"  "chr10_Persist_min2Mb.RData"
#
#


#
#
# ## FUNCTION ####################################################################
# # This function is for calculating the cross-hybridisation score between two
# # sequences.
# # getCHS <- function(
#
# K    = 7
# seq1 = seqA
# seq2 = seqB
#
#
# kmerGfree.par <- getKmerHybridisationGs(k = K, plot = FALSE)
#
# Kmers1        <- getKmers(seq.string = DNAStringSet(seq1), k = K,
# method = "Biostrings")
# RevCompKmers1 <- getKmers(seq.string=reverseComplement(DNAStringSet(seq1)), k=K,
# method = "Biostrings")
#
# Kmers2        <- getKmers(seq.string = DNAStringSet(seq2), k = K,
# method = "Biostrings")
# RevCompKmers2 <- getKmers(seq.string=reverseComplement(DNAStringSet(seq2)), k=K,
# method = "Biostrings")
# ################################################################################
#
# library(doMC)
# library(foreach)
# library(itertools)
# registerDoMC(cores=nCPU)
#
# #### FOREACH EXECUTION #########
# QP <- foreach(j=isplitVector(1:seq$length, chunks=ceiling(seq$length/SeqPartitionBy)),
# .combine="rbind", .inorder=TRUE) %dopar% {
#
#     Kmers1 <-
#
# }
# ### END OF FOREACH EXECUTION ###
#
#
#
# ################################################################################
#
#
#
#
#
#
# # Overlap with repeat Masker
# #getRepeatScores <- function(
#
#
#
#
#
#
#
#
#
# ## FUNCTION ####################################################################
#
# ################################################################################
#
# #-----------------------------------------------------------------------------
#
# ################################################################################
#
#
# ## FUNCTION ####################################################################
