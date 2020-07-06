################################################################################
# Calculate how much of the genome is covered by bins per Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
mask.dir = paste0(wk.dir, "/mask")
### OTHER SETTINGS #############################################################
chr.v <- c("chr1", "chr15")
# Permutation
NTIMES = 10000
chrlen.v = c(249250621, 102531392)
names(chrlen.v) <- chr.v
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(regioneR)
library(IRanges)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Genome
genome <- filterChromosomes(getGenome("hg19"),keep.chr=chr.v)

# Mask
mask.bed <- list.files(path=mask.dir, full.names=TRUE)
mask.bed <- sapply(X=mask.bed, simplify=FALSE, FUN=function(path){
  bed <- read.table(file=path, stringsAsFactors=FALSE, 
                    header=FALSE, row.names=NULL)[,1:3]
})
mask.bed <- do.call("rbind", mask.bed)
rownames(mask.bed) <- NULL
mask.bed <- mask.bed[mask.bed[,1]%in%chr.v,]

for(i in 1:2){
  bed <- read.table(paste0(wk.dir, "/sample_bed/samplebed", i), stringsAsFactors=FALSE)[,1:3]
  colnames(bed) <- c("chr", "start", "end")
  bed <- as.data.frame(reduce(GRanges(bed)))
  # Randomisation
  bed.r <- circularRandomizeRegions(A=bed, per.chromosome=TRUE, max.mask.overlap=1, 
                                    non.overlapping=TRUE)
  bed.r <- data.frame(bed.r, stringsAsFactors=FALSE)
  
  # Check for overlapping ranges in the randomised bed file using reduce();
  # if there no overlaps, the reduce function should have the same number of
  # ranges as the original one
  #if( nrow(bed)!=nrow(bed.r) ){ "Whole bed: Overlapping ranges in ranndomised bed." }
  #---------------------------------------
  if( !identical(nrow(bed), nrow(bed.r)) ){
    stop("Number of regions different.")
  }
  if( !identical(as.character(bed[,1]), as.character(bed.r[,1]) ) ){
    stop("Number of regions different.")
  }
  if( !identical(as.numeric(bed[,4]), as.numeric(bed.r[,4])) ){
    stop("Length distribution of regions different.") 
  }
  
  # Check inter-range length distribution of regions 
  for(chr in chr.v){
    bed.s <- bed[bed[,1]==chr,]
    i.len <- bed.s[-1,2]-bed.s[-(nrow(bed.s)),3]
    bed.r.s <- bed.r[bed.r[,1]==chr,]
    i.len.r <- bed.r.s[-1,2]-bed.r.s[-(nrow(bed.s)),3]
    #-------------------
    # Check for overlapping ranges in the randomised bed file using reduce();
    # if there no overlaps, the reduce function should have the same number of
    # ranges as the original one
    if( length( reduce(GRanges(bed.s[,1:3])) )!=length( reduce(GRanges(bed.r.s[,1:3])) ) ){
      stop("Overlapping ranges in ranndomised bed.")
    }
    #-------------------
    ind.v <- which(i.len!=i.len.r)
    for(ind in ind.v){
      if(! bed.s[ind+1,2]-bed.s[ind,3]==chrlen.v[chr]-bed.r.s[ind,3]+bed.r.s[ind+1,2] ){
        stop(paste0(chr, ": Inter-range length distribution different."))
      }
    }
    # Mismatches in calculated inter-range length distribution will arise
    # when the program reached the end of the chromosome and it has to count
    # from the start of chromosome. The inter-range length is still
    # conserved. This is why it's called 'CIRCULAR'RandomizeRegions.
    
    # Example, mismatch of inter-range length between range 450 and 451
    #> bed.s[449:451,]
    #V1        V2        V3    V4 V5
    #449 chr1 120680001 120720000 40000  *
    #450 chr1 121320001 121360000 40000  *
    #451 chr1 144000001 144040000 40000  *
    #i.len = 144000001 - 121360000 = 22640001
      
    #> bed.r.s[449:451,]
    #seqnames     start       end width strand    V4 V5
    #449     chr1 230174256 230214255 40000      * 40000  *
    #450     chr1 230814256 230854255 40000      * 40000  *
    #451     chr1   4243635   4283634 40000      * 40000  *
    #i.len.r = 249250621(*end of chr*)-230854255 = 18396366
    #i.len.r = 18396366 + 4243635 = 22640001
    
    ### i.len==i.len.r
  }
  print(i, quote=FALSE)
}

# rm(list=ls())
#---------------------------------------
# meanDistance = mean of the linear distance from every region in A to the 
# closest region in B. It is useful to answer question of the type “Are my highly 
# expressed genes closer to a certain TFBS than expected by chance?”

# Calculate distance between the two closest ends of two regions being compared
# Mean distance is always 0 if regions are next to each other or overlapping 
# (at least 1 base)
# Measures proximity to each other of regions
A <- data.frame("chr1", 11, 12)
B <- data.frame("chr1", 1, 10)
meanDistance(A,B) = 0 # next to each other
A <- data.frame("chr1", 4, 7)
B <- data.frame("chr1", 1, 10)
meanDistance(A,B) = 0 # overlapping
A <- data.frame("chr1", 13, 16)
B <- data.frame("chr1", 1, 10)
meanDistance(A,B) = 2 # length of intervening bases between end of B and start of A

#---------------------------------------
# numOverlaps = number of regions in A overlapping at least one region in B
# Does not reduce any of A or B so I can represent the actual number of bins in
# the B bed file
genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
A <- toGRanges( data.frame("chr1", c(1, 10, 20, 27, 26), c(12, 13, 28, 40, 39)) )
B <- toGRanges(data.frame("chr1", c(25, 27, 26), c(34, 35, 35)))
A.red <- reduce( toGRanges( data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40)) ) )
B.red <- reduce( toGRanges(data.frame("chr1", c(25, 27), c(35, 35))) )
# A - B = 3
numOverlaps(A, B, count.once=TRUE)
# B - A = 3
numOverlaps(B, A, count.once=TRUE)
# A.red - B = 2
numOverlaps(A.red, B, count.once=TRUE)
# B.red - A = 1
numOverlaps(B.red, A, count.once=TRUE)
# A.red - B.red = 2
numOverlaps(A.red, B.red, count.once=TRUE)
# B.red - A.red = 1
numOverlaps(B.red, A.red, count.once=TRUE)
#---------------------------------------
# Custom evaluation function = calculate overlap percentage relative to total length of regions
# to differentiate degrees of overlap (amount of bases overlapping)

A <- data.frame("chr1", c(1, 10, 20, 30), c(12, 13, 28, 40))
B <- data.frame("chr1", 25, 35)
percOlap <- function(A, B, ...) {
  total <- sum( c(width(A), width(B)) )
  common <- sum( width(commonRegions(A,B)) )
  return( common/total*100 )
}
commonRegions(A,B)
numOverlaps(A, B)
numOverlaps(A, B, count.once=TRUE)

A <- data.frame("chr1", 12, 12)
B <- data.frame("chr1", 1, 10)
meanDistance(A,B)
A <- data.frame("chr1", 4, 7)
B <- data.frame("chr1", 1, 10)
meanDistance(A,B)
#---------------------------------------
A <- toGRanges( data.frame("chr1", c(1, 5, 20, 30), c(8, 13, 28, 40), x=c(1,2,3,4),
                y=c("a", "b", "c", "d")) )
B <- toGRanges( data.frame("chr1", 25, 35) )
bed <- overlapRegions(A, B, type="any")[,c("chr", "startA", "endA")]
#---------------------------------------
# Local z-score
genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
A <- createRandomRegions(nregions=20, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))

pt <- overlapPermTest(A=A, B=B, ntimes=10, genome=genome, non.overlapping=FALSE)
plot(pt)

lz <- localZScore(A=A, B=B, pt=pt)
plot(lz)


pt2 <- permTest(A=A, B=B, ntimes=10, randomize.function=randomizeRegions, evaluate.function=list(overlap=numOverlaps, distance=meanDistance), genome=genome, non.overlapping=FALSE)
plot(pt2)

lz2 <- localZScore(A=A, B=B, pt2)
plot(lz2)

# rm(list=ls())
