## FUNCTION ####################################################################
HicHybridisation <- function(
  # Number of CPU cores to be used
  nCPU = 2, # optimal for the task on bucephalus (takes ~5GB per core for chr1)
  # Chromosome for which the data are to be pulled:
  chr = "chr1",
  # Feature database path and filename suffix:
  featureDBPath = #"/home/alex/Desktop/CHROMSEQ/OUT",
        "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc",
  # Feature database filename suffix
  suffix = "min2Mb",
  k = 7,
  scale = 40000 # the length of DNA/Hi-C resolution for Gfree scaling to per-nt
){
################################################################################
library(Biostrings)
library(RColorBrewer)
library(doMC)
library(foreach)
library(itertools)
registerDoMC(cores=nCPU)
################################################################################

id <- paste0(chr,"_",suffix,"_kmer",k)

BINKMER.name <- load(paste0(featureDBPath,"/",chr,"_BinKmer",k,"_",suffix,".RData"))
eval(parse(
 text=paste0("BINKMER.MX <- ",BINKMER.name,"; rm(",BINKMER.name,"); rm(BINKMER.name); gc()")
))
load(paste0(featureDBPath,"/",chr,"_Persist_",suffix,".RData"))
gfree <- read.table(file=paste0("./Gfree_",k,"mer.par"), header=TRUE, stringsAsFactors=FALSE)

############################
kmer.vec <- dimnames(BINKMER.MX)[[2]][-c(1:3)]
kmer.vec.revcomp <- as.vector(reverseComplement(DNAStringSet(kmer.vec)))
kmer.revcomp.ind <- match(kmer.vec, kmer.vec.revcomp)
rm(kmer.vec, kmer.vec.revcomp)
############################
hits.len <- length(PERSIST.MX$hits[,1])
ij.mx <- cbind(i=PERSIST.MX$hits[,"i"], j=PERSIST.MX$hits[,"j"])
ntis.vec <- PERSIST.MX$ntis
foreach.export <- c("BINKMER.MX","ij.mx","kmer.revcomp.ind","gfree","scale")
rm(PERSIST.MX)
############################

gc()

#### FOREACH EXECUTION #########
HYB.MX <- foreach(itr=isplitVector(1:hits.len, chunks=nCPU),
                      .combine="cbind", .inorder=TRUE,
                      .export=foreach.export,
                      .noexport=ls()[!ls()%in%foreach.export]
                     ) %dopar% {

  hyb.data <- sapply(itr,
   FUN=function(i){
     kmer.mx <- BINKMER.MX[which(BINKMER.MX[,"bins"] %in% ij.mx[i,]), -c(1:3) ]
     #kmerA <- kmer.mx[1,]
     #kmerB <- kmer.mx[2,]
     #kmerArc <- kmer.mx[1,kmer.revcomp.ind]
     #kmerBrc <- kmer.mx[2,kmer.revcomp.ind]
     kmerAall <- (kmer.mx[1,] + kmer.mx[1,kmer.revcomp.ind]) #(kmerA + kmerArc)
     kmerBall <- (kmer.mx[2,] + kmer.mx[2,kmer.revcomp.ind]) #(kmerB + kmerBrc)
     return( c(
       # everything is counted twice, hence 2 in the division
       sum(pmin(kmerAall, kmerBall)*gfree[,2])/(2*scale),
       # everything is twice in frequency
       #hist(kmerAall - kmerBall, breaks=100, col="navy")
       sd( kmerAall - kmerBall )
       #mean( kmerAall - kmerBall ) # mean is alws 0; double checked numerically
     ) )
   }, simplify=TRUE) # USE.NAMES=FALSE

   return(hyb.data)

}
### END OF FOREACH EXECUTION ###


rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))

df <- data.frame(ntis=ntis.vec, Gfree=HYB.MX[1,])
pdf(file=paste0(id,"_Gfree_boxplot.pdf"), width=9, height=7)
  boxplot(formula=Gfree~ntis, data=df,
          col=rf(max( ntis.vec )), cex=0.1, # ylim=c(0,30)
          outline=FALSE, xlab="Num(non-0 tissues/cell lines)",
          ylab="T.H.E. Gfree (kcal/mol)",
          main=paste0(id))
dev.off()

df <- data.frame(ntis=ntis.vec, SD=HYB.MX[2,])
pdf(file=paste0(id,"_Discordance_boxplot.pdf"), width=9, height=7)
  boxplot(formula=SD~ntis, data=df,
          col=rf(max( ntis.vec )), cex=0.1, # ylim=c(0,30)
          outline=FALSE, xlab="Num(non-0 tissues/cell lines)",
          ylab="T.H.E. sd(discordance)",
          main=paste0(id))
dev.off()

############################
dimnames(HYB.MX)[[1]] <- c("Gfree", "Discordance")
save(HYB.MX, file=paste0(featureDBPath,"/",chr,"_Hyb",k,"_",suffix,".RData"))
############################

rm(rf, df, HYB.MX); gc()
print("HicHybridisation is DONE!", quote=FALSE)

}
################################################################################
suppressWarnings(suppressPackageStartupMessages(library(compiler)))
HicHybridisation <- cmpfun(HicHybridisation, options=list(suppressUndefined=TRUE))

# Execution of the function written above...
################################################################################
HicHybridisation(chr="chr1",suffix="min2Mb",k=7,nCPU=4,scale=40000,
                 featureDBPath = "/home/alex/Desktop/CHROMSEQ/OUT")

HicHybridisation(chr="chr1",suffix="min05Mb",k=7,nCPU=4,scale=40000,
                 featureDBPath = "/home/alex/Desktop/CHROMSEQ/OUT")

HicHybridisation(chr="chr1",suffix="min2Mb",k=4,nCPU=4,scale=40000,
                 featureDBPath = "/home/alex/Desktop/CHROMSEQ/OUT")

HicHybridisation(chr="chr1",suffix="min05Mb",k=4,nCPU=4,scale=40000,
                 featureDBPath = "/home/alex/Desktop/CHROMSEQ/OUT")
