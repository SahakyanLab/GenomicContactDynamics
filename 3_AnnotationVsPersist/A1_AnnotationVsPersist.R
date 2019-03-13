# AnnotationVsPersist
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
                            # "AlexMac", "AlexCluster"

#*** See whether the content of the <lib>, specified below, can be moved to the
#*** actual project folder for self consistency. Then, specify the location of
#*** <lib> within the whorunsit block below. If you put the dependent function
#*** in <lib> into the <objective.dir>/lib location, then the following line
#*** outside the whorunsit block, after output.dir specification is enough:
#***    lib = paste0(objective.dir, "/lib")

# Functions to source
# lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelCluster"){
    # Main directory for the task
    objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
    # Annotation file
    annoFile.dir   = "/t1-home/icbml/ltamon/Database/ucsc_tables"
    # HiC_Human21 persistent contact files
    persist.dir   = "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
    os = "Linux"
  } else if(whorunsit == "LiezelMac"){
    objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
    annoFile.dir   = "/Users/ltamon/Database/ucsc_tables"
    persist.dir   = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
    os = "Mac"
  } else if(whorunsit == "AlexMac"){
    objective.dir = "./"
    annoFile.dir   = "/Volumes/Data/Database/ucsc_tables"
    persist.dir   = "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc"
    os = "Mac"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}

output.dir = paste0(objective.dir, "/out_AnnotationVsPersist")
### OTHER SETTINGS #############################################################
# "2"(2MB gap) or "05"(0.5 MB minimum gap), refers to a minimum gap accepted to
# further process a given contact: two points should be far enough to filter for
# contacts within a TAD.
gc.v = c("2")
chr.v = c(1:22, "X")

# Length of ij bins
bin.length = 40000L

# Genome version
genome.ver = "hg19"
anno.filename = paste0(genome.ver, "anno")

# Set overlap minimum and gap maximum for overlap of transcript/annotations
# with ij bins
olap.min = 1L
#*** os is defined in whorunsit block
if( os == "Linux" ){
  gap.max = 0L
} else if( os == "Mac" ){
  gap.max = -1L
}
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(data.table)
library(IRanges)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################

# load annotation file
anno.file <- fread(file=paste0(annoFile.dir, "/", anno.filename),
                   header=TRUE, data.table=FALSE)
# check for missing values in annotation file
if( sum(apply(anno.file, MARGIN=c(1,2), function(x) any(is.na(x))))!=0 ){
  stop("Missing values in annotation file.")
}
# lengths of transcripts
txSize <- anno.file[,"txEnd"]-anno.file[,"txStart"]
#just checking if txEnd always greater than txStart as assumed
#in the annotation file
if( sum(txSize<0) ){
  stop("There are txStart > txEnd.")
}

for(gc in gc.v){
  # work per chromosome, subset annotation file and load persist matrix per chromosome
  ANNO.PERSIST.MX <- foreach(chr=chr.v, .inorder=TRUE, .combine="rbind") %do% {
    load(file=paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))
    # get the relevant bins (from persist matrix) for that chromosome
    # bins sorted to increasing values so the ranges are also increasing
    bins.uniq <- unique( c(unique(PERSIST.MX$hits[,"i"]),
                           unique(PERSIST.MX$hits[,"j"])) )
    bins.uniq <- sort(x=bins.uniq, decreasing=FALSE)
    bin.end   <- bins.uniq*bin.length
    # subset annotation file per chromosome
    subset.anno <- subset(x=anno.file, chrom==paste0("chr", chr))
    # assign bins based on txStart-txEnd
    # query   <- subset.anno; annotation (to be classified)
    # subject <- bins.uniq (sorted, increasing); bins (classification)
    olap.mx <- WhichOverlap(start.query=subset.anno[,"txStart"],
                            end.query=subset.anno[,"txEnd"],
                            space.query=rep("a", nrow(subset.anno)),
                            start.subject=bin.end-bin.length+1,
                            end.subject=bin.end,
                            space.subject=rep("a",length(bin.end)),
                            maxgap=gap.max, minoverlap=olap.min)

    # annotation assigned to bin
    # note that one annotation can be assigned to multiple bins
    col12 <- cbind(subset.anno[olap.mx[,"query"],],
                   HiC21bin=bins.uniq[olap.mx[,"subject"]])
    # select all contact pairs containing the annotation (either in i or j bin)
    # make sure indices are unique
    lst <- list()
    for( bin in unique(col12[, "HiC21bin"]) ){
      # assumes that the contacts from PERSIST.MX are unique
      jbinpartnerOfBin <- PERSIST.MX$hits[PERSIST.MX$hits[,"i"]==bin,"j"]
      ibinpartnerOfBin <- PERSIST.MX$hits[PERSIST.MX$hits[,"j"]==bin,"i"]
      # all the indices that are paired with bin
      hits.noNA.ij <- c(jbinpartnerOfBin, ibinpartnerOfBin)
      #indices of ij contacts in persist matrix that has the bin of interest
      ind.all <- unique(c(which(PERSIST.MX$hits[,"i"]==bin),
                          which(PERSIST.MX$hits[,"j"]==bin)))
      #hits.NA      <- apply(persist.mx[ind.all,c(i, j)], MARGIN=c(1,2),
      #                      function(x) ifelse(x==bin, x<-NA, x<-x))
      #hits.noNA.ij <- na.omit(c(rbind(hits.NA[,i], hits.NA[,j])))
      strengthString <- apply(X=PERSIST.MX$hits[ind.all,3:23], MARGIN=1,
                              FUN=function(x){paste(x, collapse=";")})
      tissuenames <- paste(colnames(PERSIST.MX$hits[ind.all,3:23]),
                           collapse="_")
      all.val <- cbind(HiC21binPartner=hits.noNA.ij,
                       ntis=PERSIST.MX$ntis[ind.all],
                       strengthString,
                       SumStrength=PERSIST.MX$valsum[ind.all],
                       SumStrengthOverNtis=PERSIST.MX$valsum[ind.all]/PERSIST.MX$ntis[ind.all])
      setnames(as.data.frame(all.val), old=3, new=tissuenames)
      lst[[as.character(bin)]] <- c( HiC21bin=bin,
                                     apply(all.val, MARGIN=2, function(x){paste(x, collapse=",")}) )
    }
    ANNO.PERSIST.MX <- merge(col12, do.call(rbind, lst),
                             by.x="HiC21bin", all.x=TRUE) #keep rows of x
    save(ANNO.PERSIST.MX, file=paste0(output.dir, "/", "chr", chr, "_min", gc, "Mb_",
                                      genome.ver, "_ALL_AnnoVsPersist_olapMin",
                                      olap.min, ".RData"))
    ANNO.PERSIST.MX
  }
  save(ANNO.PERSIST.MX, file=paste0(output.dir, "/", "chrALL_min", gc, "Mb_",
                                    genome.ver, "_ALL_AnnoVsPersist_olapMin",
                                    olap.min, ".RData"))
}

#rm(list=ls())
