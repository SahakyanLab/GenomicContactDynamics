################################################################################
# Using the (1) control matrix with all the HiC bins and (2) previously-generated 
# UNIQBIN.DF that contains data on unique bins per Cp per chr, obtain object KMERCP, 
# a list of 
# a. mean matrix (7mer X Cp) - mean of fractions per k-mer per Cp, NAs removed
# during calculation
# b. sd matrix (7mer X Cp) - sd of fractions per k-mer per Cp, NAs removed
# during calculation
# c. neglopg10pval matrix - from Mann–Whitney–Wilcoxon (MWW) (non-parametric, two 
# independent distributions) to compare, for each k-mer, its fraction values per Cp
# to that of the control. Non-parametric to better account for difference in
# number of bins contributing to the test and control set. 

# The MWW test is directional. If mean.test > mean.control, alternative hypothesis
# is greater. If mean.test < mean.control, alternative hypothesis is less.
# If mean.test = mean.control, alternative hypothesis is two.sided.

# When p-value=0, -log10(p-value)=+Inf. For plotting, this is converted to
# the maximum value of -log10(p-value).  
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    # UNIQBIN.DF directory
    data.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/9_GenomicComposition"
  } else if(whorunsit == "LiezelCluster"){
    # UNIQBIN.DF directory
    data.dir = "/t1-data/user/ltamon/Database/HiC_features_GSE87112_RAWpc"
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/9_GenomicComposition"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
control.dir = paste0(wk.dir, "/out_control")
out.dir = paste0(wk.dir, "/out_compare")
### OTHER SETTINGS #############################################################
gcb = "min05Mb"
kmer.len = "7"
# Number of kmers
nCPU = 1L #~100G
kmerDistVal = "Median"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cp.v <- as.character(1:21)

# Load KMERfrBIN.HiCAll (control)
load(file=paste0(control.dir, "/HiCAll_KMERfrBINPCP", kmer.len, "_", gcb, ".RData"))
kmer.v <- dimnames(KMERfrBIN.HiCAll)[[1]]
kmer.v.len <- length(kmer.v)

# Initialize output matrices
KMERCP <- list()
# Plus 1 column for the HiC all control
KMERCP$CENTR <- matrix(data=NA, nrow=kmer.v.len, ncol=length(cp.v)+1L)
dimnames(KMERCP$CENTR) <- list(kmer.v, c("HiCAll", cp.v))
KMERCP$SD <- KMERCP$CENTR
KMERCP$NEGLOG10PVAL <- KMERCP$CENTR[,-1]

# Add control values
if(kmerDistVal=="Mean"){
  
  KMERCP$CENTR[,"HiCAll"] <- rowMeans(KMERfrBIN.HiCAll, na.rm=TRUE)
  print(paste0("Central value is MEAN."), quote=FALSE)
  
} else if(kmerDistVal=="Median"){
  
  KMERCP$CENTR[,"HiCAll"] <- apply(X=KMERfrBIN.HiCAll, MARGIN=1, FUN=function(x){
    median(x=x, na.rm=TRUE)
  })
  print(paste0("Central value is MEDIAN."), quote=FALSE)
  
} else {
  
  stop("Checkpoint 1.")
  
}

KMERCP$SD[,"HiCAll"] <- apply(X=KMERfrBIN.HiCAll, MARGIN=1, FUN=function(x){
  sd(x=x, na.rm=TRUE)
})

# Load UNIQBIN.DF
load(file=paste0(data.dir, "/chrALL_Uniqbin_", gcb, ".RData"))

UNIQBIN.DF <- by(data=paste(UNIQBIN.DF[,"chr"], UNIQBIN.DF[,"bin"], 
                            sep="_"),
                 INDICES=UNIQBIN.DF[,"cp"],
                 FUN=function(x){
                   as.character(x)
                 })

for(cp in cp.v){
  
  print(paste0("cp=", cp, "..."), quote=FALSE)
  
  KFB.log <- colnames(KMERfrBIN.HiCAll)%in%UNIQBIN.DF[[cp]]
  
  toExport <- c("kmer.v", "KFB.log", "KMERfrBIN.HiCAll")
  
  #### PARALLEL EXECUTION #########
  
  MX <- foreach(itr=isplitVector(1:kmer.v.len, chunks=nCPU), 
                .inorder=TRUE, .combine="rbind",
                .export=toExport, .noexport=ls()[!ls()%in%toExport]
                
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      kmer <- kmer.v[i]
      
      if(kmerDistVal=="Mean"){
        samp.x.centr <- mean(KMERfrBIN.HiCAll[kmer,KFB.log], na.rm=TRUE)
        cont.y.centr <- mean(KMERfrBIN.HiCAll[kmer,], na.rm=TRUE)
      } else if(kmerDistVal=="Median"){
        samp.x.centr <- median(KMERfrBIN.HiCAll[kmer,KFB.log], na.rm=TRUE)
        cont.y.centr <- median(KMERfrBIN.HiCAll[kmer,], na.rm=TRUE)
      } else {
        stop("Checkpoint 2.")
      }
      
      return.v <- c(samp.x.centr, 
                    sd(KMERfrBIN.HiCAll[kmer,KFB.log], na.rm=TRUE))
      
      if(samp.x.centr==cont.y.centr){
        
        mw <- wilcox.test(x=as.vector( na.omit(KMERfrBIN.HiCAll[kmer,KFB.log]) ), 
                          y=as.vector( na.omit(KMERfrBIN.HiCAll[kmer,]) ) ,
                          alternative="two.sided", paired=FALSE)
        
      } else {
        
        if(samp.x.centr>cont.y.centr){
          alt.hyp <- "greater"
        } else if(samp.x.centr<cont.y.centr){
          alt.hyp <- "less"
        } else {
          stop("Checkpoint 1.")
        } 
        
        # Wilcoxon-Mann-Whitney test (non-parametric, unpaired/independent groups) 
        # One-sided alternative "greater" is that x is shifted to the right of y,
        # alternative=greater, x>y?
        mw <- wilcox.test(x=as.vector( na.omit(KMERfrBIN.HiCAll[kmer,KFB.log]) ), 
                          y=as.vector( na.omit(KMERfrBIN.HiCAll[kmer,]) ) ,
                          alternative=alt.hyp, paired=FALSE)
      }
      
      # If p.value=0, log10 would give positive Inf
      return( c(return.v, -(log10(mw$p.value))) ) 
      
    }) # sapply chunk end
    
    return(do.call("rbind", chunk))
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  ind <- which(colnames(KMERCP$CENTR)==cp)
  KMERCP$CENTR[,ind] <- MX[,1]
  KMERCP$SD[,ind] <- MX[,2]
  # Because it has one less column (no HiCAll)
  KMERCP$NEGLOG10PVAL[,ind-1L] <- MX[,3]
  
  rm(MX, KFB.log); gc()
  
  print(paste0("cp=", cp, " done!"), quote=FALSE)
  
} # cp.v for loop end

save(KMERCP, file=paste0(out.dir, "/KMERCP", kmer.len, "_", gcb, "_", 
                         kmerDistVal, ".RData"))

# rm(list=ls()) 