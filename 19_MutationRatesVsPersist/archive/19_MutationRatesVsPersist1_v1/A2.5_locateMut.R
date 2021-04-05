################################################################################
# After liftover, check if the hg19 non-coding mutations overlap with 
# (a) coding transcript exon UTRs, (b) coding transcript translated exons 
# (c) non-coding transcript exons. Also, check all mutations are within hg19 
# chromosome length. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
mut.dir = paste0(wk.dir, "/out_filter")
out.dir = paste0(wk.dir, "/out_locateMut")
bed.dir = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/regions_anno/out_extractRegions")
chrLenFile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
LO.chain = "Hg19" # "Hg19" | "hg38ToHg19"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(GenomicRanges)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
#fileNme <- paste0("CosmicNCV_", LO.chain, "_final_2000.RData") # REMOVE, ONLY FOR TESTING
fileNme <- paste0("CosmicNCV_", LO.chain, "_final.RData")
load(file=paste0(mut.dir, "/", fileNme))

# Output of liftover can have missing coordinates. These probably yielded two or more liftOver 
# coordinates per mutation.
if( any(is.na(ncv.df$start)) ){
  print("Missing start.")
}
if( any(is.na(ncv.df$end)) ){
  print("Missing end.")
}

temp <- c(ncv.df$start, ncv.df$end)
if( any( !is.na(temp) & temp<=0 )  ){
  stop("Negative or zero start and end coordinates")
}

# Check if mutations are within hg19 chromosome length
chrLen.df <- read.delim(file=chrLenFile, header=TRUE, stringsAsFactors=FALSE)
chr.v <- unique(ncv.df$chr)
outChr.TF <- unlist(sapply(X=na.omit(chr.v), simplify=FALSE, FUN=function(chr){
  chr.len <- chrLen.df[chrLen.df$chromosome==chr,"length.bp"]
  return( max(ncv.df[ncv.df$chr==chr,"end"], na.rm=TRUE)>chr.len )
}))
print(outChr.TF)
outChr.TF <- unname(na.omit(outChr.TF))
if( any(outChr.TF) ){
  stop("Mutations outside hg19 chromosome length.")
}
rm(chrLen.df, outChr.TF)

# Check overlap with beds of interests

for( bedNme in list.files(bed.dir, full.names=FALSE) ){
  
  bed.gr <- read.table(file=paste0(bed.dir, "/", bedNme), header=FALSE, sep="\t", 
                       stringsAsFactors=FALSE)
  bed.gr <- GRanges(IRanges(start=bed.gr$V2, end=bed.gr$V3), seqnames=bed.gr$V1, 
                    strand="*", score=bed.gr$V4) 
  mut.gr <- GRanges(IRanges(start=ncv.df$start, end=ncv.df$end), seqnames=ncv.df$chr, 
                    strand="*", score=1:nrow(ncv.df)) 
  countOlap <- GenomicRanges::countOverlaps(query=mut.gr, subject=bed.gr, minoverlap=1L, 
                                            maxgap=-1L, type="any")
  Olap <- GenomicRanges::findOverlaps(query=mut.gr, subject=bed.gr, minoverlap=1L, 
                                            maxgap=-1L, type="any")
  
  if( length(countOlap)!=length(mut.gr) ){
    stop("Checkpoint 1.")
  }
  
  temp <- table(queryHits(Olap))
  if( any(countOlap[as.numeric(names(temp))]!=unname(temp)) ){
    stop("Checkpoint 2.")
  }
  rm(temp, Olap)
  
  countOlap <- sum(countOlap>0)
  
  if(countOlap>0){
    write(x=paste0(countOlap, " mutations (out of ", nrow(ncv.df), ") in ", 
                   fileNme, " overlap with ", bedNme), append=TRUE,
          file=paste0(out.dir, "/", LO.chain, "_output.txt"))
  }
  
  rm(bed.gr, mut.gr, countOlap)
  
}

# rm(list=ls()); gc()
