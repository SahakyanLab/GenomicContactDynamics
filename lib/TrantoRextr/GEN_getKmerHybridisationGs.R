## FUNCTION ####################################################################
# For complete DNA k-meric set of a given k, where k>2, generates a table of   #
# all sequence permutations and the estimates of full hybridisation energies   #
# based on the triplet hybridisation energies (by sliding window + averaging). #
# REQUIRES:                                                                    #
# library(gtools)
# library(Biostrings)
# LIB.TRANTOR = "/Users/alex/GIT/GITrepo/TrantoR"
# source(paste0(LIB.TRANTOR, "/GEN_getKmers.R"))
################################################################################
getKmerHybridisationGs <- function(k = 7, plot = FALSE){

  triad.parameters <- c("AAA,-0.844597",
                        "AAC,-1.841904",
                        "AAG,-1.201194",
                        "AAT,-0.991596",
                        "ACA,-1.121939",
                        "ACC,-1.793995",
                        "ACG,-1.615048",
                        "ACT,-0.781693",
                        "AGA,-1.103536",
                        "AGC,-1.528461",
                        "AGG,-1.323278",
                        "ATA,-0.562379",
                        "ATC,-1.157521",
                        "ATG,-1.263601",
                        "CAA,-0.988509",
                        "CAC,-2.088824",
                        "CAG,-1.625284",
                        "CCA,-1.568813",
                        "CCC,-2.396507",
                        "CCG,-1.888906",
                        "CGA,-1.668273",
                        "CGC,-2.195726",
                        "CTA,-0.871636",
                        "CTC,-1.19845",
                        "GAA,-1.317278",
                        "GAC,-1.498999",
                        "GCA,-1.45443",
                        "GCC,-1.973081",
                        "GGA,-1.696158",
                        "GTA,-1.158422",
                        "TAA,-0.519499",
                        "TCA,-1.042342")
  triad.parameters <- read.csv(textConnection(triad.parameters),
                               header=FALSE, as.is=TRUE)

  #-----------------------------------------------------------------------------
  library(gtools)
  #-----------------------------------------------------------------------------
  query.permuts.tri <- permutations(n=4, r=3, v=c("A","G","T","C"),
                                    repeats.allowed=TRUE)
  #-----------------------------------------------------------------------------
  query.permuts.tri <- sapply(1:length(query.permuts.tri[,1]), FUN=function(i){
    return(paste(query.permuts.tri[i,],collapse=""))
  }, simplify=TRUE, USE.NAMES=FALSE)

  sorted.Gfrees.triad <- NULL
  for(j in query.permuts.tri){
    ind <- which(triad.parameters[,1]==j)
    if(length(ind)==0){
      ind <- which(triad.parameters[,1]==
                   as.character(reverseComplement(DNAString(j))))
    }
    sorted.Gfrees.triad <- c(sorted.Gfrees.triad, triad.parameters[ind,2])
  }
  # paste(query.permuts.tri, sorted.Gfrees.triad)
  # sorted.Gfrees.triad - Gfree parameters for triads retrieved for all 64
  # permutations (out of the 32 unique ones in triad.parameters), sorted in
  # the same order (query.permuts.tri) as the 3-meric outout of getKmers()
  # with k = 3.
  #-----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  query.permuts <- permutations(n=4, r=k, v=c("A","G","T","C"),
                                repeats.allowed=TRUE)
  #-----------------------------------------------------------------------------
  query.permuts <- sapply(1:length(query.permuts[,1]), FUN=function(i){
    return(paste(query.permuts[i,],collapse=""))
  }, simplify=TRUE, USE.NAMES=FALSE)

  Gfrees <- rep(NA, length(query.permuts))

  for(i in 1:length(query.permuts)){
    #print(i)
    kmer <- query.permuts[i]
    triads <- getKmers(seq.string=DNAStringSet(kmer), k=3, method="Biostrings")
    Gfrees[i] <- sum(triads*sorted.Gfrees.triad)/(k-2)
  }

  if(plot){
    pdf(height=5, width=8, file=paste0("hist_AllUnique_k",k,".pdf"))
        hist(Gfrees, breaks=50, col="cornflowerblue", xlim=c(-2.5,-0.5),
             xlab="dG score, kcal/mol", ylab="count", main="")
    dev.off()
  }
  #-----------------------------------------------------------------------------

  print(paste0("getKmerHybridisationGs: returning Gfree table for all ",
               k,"-mers."), quote=FALSE)
  return(
    data.frame(kmers=query.permuts, Gfrees=Gfrees, stringsAsFactors=FALSE)
  )
}
################################################################################
