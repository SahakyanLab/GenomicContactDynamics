## FUNCTION ####################################################################
# This function takes the chromosome identifier, Hi-C experiment resolution (i.e.
# the real length of its bins), the vector of unique bin IDs for which to examine
# the overlaps and the Repeat Masker database location. The program returns a
# matrix with bin number, startpos, endpos and the counts of all overlaped repeat
# families/subfamilies.
################################################################################
getBinsOverlapRepMasker <- function(
  chr = "chr1",
  hic.resol = 40000,
  # unique Hi-C bin ids
  bins = unique( c(unique(PERSIST.MX$hits[,"i"]),
                   unique(PERSIST.MX$hits[,"j"])) ),
  repmask.Filepath = paste0(data.dir, "/ucsc_tables/hsa_RepeatMasker/RepeatMasker_hg19"),
  overlapBy = "repName",
  min.olap = min.olap,
  max.gap = max.gap,
  type.olap = type.olap
){
  repmask <- getRepMaskSubset(chr=chr, freadchar=" ", 
                              repmask.Filepath=repmask.Filepath,
                              repmask.col=c("genoName", "genoStart","genoEnd","strand",
                                            "repName","repClass","repFamily"))
  
  # Ignores 2 repNames with two families/classes designation because they're not 
  # transposons but simple repeat/satellite
  all.overlapBy <- unique(repmask[[overlapBy]])
  
  bins.len <- length(bins)
  endpos   <- bins*hic.resol
  startpos <- endpos-hic.resol+1
  
  repeat.counts <- matrix(0, nrow=bins.len, ncol=length(all.overlapBy))
  dimnames(repeat.counts)[[2]] <- all.overlapBy
  rm(all.overlapBy)
  
  overlap <- WhichOverlap(start.query   = startpos,
                          end.query     = endpos,
                          space.query   = rep(chr, length(startpos)),
                          start.subject = repmask$genoStart,
                          end.subject   = repmask$genoEnd,
                          space.subject = rep(chr, length(repmask$genoStart)),
                          minoverlap    = min.olap,
                          maxgap        = max.gap,
                          type          = type.olap)
  
  for(k in 1:bins.len){
    ovr.ind <- which(overlap[,"query"]==k)
    if(length(ovr.ind)!=0){
      reps.k <- table(repmask[unique(overlap[ovr.ind,"subject"]), overlapBy])
      repeat.counts[k,names(reps.k)] <- reps.k
    }
  }
  
  return(cbind(bins=bins, startpos=startpos, endpos=endpos, repeat.counts))
  
}
################################################################################
