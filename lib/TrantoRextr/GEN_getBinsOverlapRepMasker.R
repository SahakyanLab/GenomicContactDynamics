## FUNCTION ####################################################################
# This function takes the chromosome identifier, Hi-C experiment resolution (i.e.
# the real length of its bins), the vector of unique bin IDs for which to examine
# the overlaps and the Repeat Masker database location. The program returns a
# matrix with bin number, startpos, endpos and the counts of all overlaped repeat
# families.
################################################################################
getBinsOverlapRepMasker <- function(
  chr="chr10",
  hic.resol=40000,
  # unique Hi-C bin ids
  bins=unique(as.vector(as.matrix(PERSIST.MX$hits[,c("i","j")]))),
  repmask.Filepath="/Volumes/Data/Database/RepeatMasker_hg19/RepeatMasker_hg19"
){
  
  repmask <- getRepMaskSubset(chr=chr, freadchar="\t", 
                              repmask.Filepath=repmask.Filepath,
                              repmask.col=c("genoStart","genoEnd","strand","repName","repClass","repFamily"))
  
  all.repFamily <- unique(repmask$repFamily)
  
  bins.len <- length(bins)
  startpos <- (bins*hic.resol)-hic.resol+1
  endpos   <- bins*hic.resol
  
  repeat.counts <- matrix(0, nrow=bins.len, ncol=length(all.repFamily))
  dimnames(repeat.counts)[[2]] <- all.repFamily
  
  overlap <- WhichOverlap(start.query   = startpos,
                          end.query     = endpos,
                          space.query   = rep(chr, length(startpos)),
                          start.subject = repmask$genoStart,
                          end.subject   = repmask$genoEnd,
                          space.subject = rep(chr, length(repmask$genoStart)),
                          maxgap        = -1L)
  
  for(k in 1:bins.len){
    #print(k)
    ovr.ind <- which(overlap[,"query"]==k)
    if(length(ovr.ind)!=0){
      reps.k <- table(repmask[unique(overlap[ovr.ind,"subject"]),"repFamily"])
      repeat.counts[k,names(reps.k)] <- reps.k
    }
  }
  
  return(cbind(bins=bins, startpos=startpos, endpos=endpos, repeat.counts))
  
}
################################################################################
