################################################################################
#AnnotationVsPersist

################################################################################ 
#functions to source
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
#annotation file
#annoFile.dir = "/Users/ltamon/Database/ucsc_tables"
annoFile.dir = "/t1-home/icbml/ltamon/Database/ucsc_tables"
#HiC_Human21 persist files
#persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
persist.dir = "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
#output.dir = paste0(objective.dir, "/out_AnnotationVsPersist_test")    
output.dir = paste0(objective.dir, "/out_AnnotationVsPersist")    

library(foreach)
library(data.table) 
library(IRanges) 

#length of ij bins
bin.length = 40000L

#genome version
genome.ver = "hg19"
anno.filename = paste0(genome.ver, "anno") 

#2(2MB gap) or "05"(0.5 MB gap), refers to minimum gap between contacting bins, 
#two points should be far enough to filter for contacts within a TAD
gc.v = c("2")

chr.v = c(1:22, "X") 

#set overlap minimum and gap maximum for overlap of transcript/annotations 
#with ij bins
olap.min = 1L 
gap.max = 0L #default for MAC (-1L) for Linux (0L)

source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))

################################################################################ 
################################################################################ 
#load annotation file 
anno.file <- fread(file=paste0(annoFile.dir, "/", anno.filename), 
                   header=TRUE, data.table=FALSE)
#check for missing values in annotation file
if( sum(apply(anno.file, MARGIN=c(1,2), function(x) any(is.na(x))))!=0 ){
  stop("Missing values in annotation file.")
}
#lengths of transcripts
txSize <- anno.file[,"txEnd"]-anno.file[,"txStart"]
#just checking if txEnd always greater than txStart as assumed 
#in the annotation file
if( sum(txSize<0) ){
  stop("There are txStart > txEnd.")
} 

for(gc in gc.v){
  #work per chromosome, subset annotation file and load persist matrix per chromosome
  ANNO.PERSIST.MX <- foreach(chr=chr.v, .inorder=TRUE, .combine="rbind") %do% {
    load(file=paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))
    #get the relevant bins (from persist matrix) for that chromosome
    #bins sorted to increasing values so the ranges are also increasing
    bins.uniq <- unique( c(unique(PERSIST.MX$hits[,"i"]),
                           unique(PERSIST.MX$hits[,"j"])) )
    bins.uniq <- sort(x=bins.uniq, decreasing=FALSE)
    bin.end   <- bins.uniq*bin.length
    #subset annotation file per chromosome
    subset.anno <- subset(x=anno.file, chrom==paste0("chr", chr))
    #assign bins based on txStart-txEnd
    #query   <- subset.anno; annotation (to be classified)
    #subject <- bins.uniq (sorted, increasing); bins (classification)
    olap.mx <- WhichOverlap(start.query=subset.anno[,"txStart"], 
                            end.query=subset.anno[,"txEnd"], 
                            space.query=rep("a", nrow(subset.anno)),
                            start.subject=bin.end-bin.length+1, 
                            end.subject=bin.end, 
                            space.subject=rep("a",length(bin.end)),
                            maxgap=gap.max, minoverlap=olap.min)
    #annotation assigned to bin
    #note that one annotation can be assigned to multiple bins
    col12 <- cbind(subset.anno[olap.mx[,"query"],], 
                   HiC21bin=bins.uniq[olap.mx[,"subject"]])
    #select all contact pairs containing the annotation (either in i or j bin)
    #make sure indices are unique 
    lst <- list()
    for( bin in unique(col12[, "HiC21bin"]) ){
      #assumes that the contacts from PERSIST.MX are unique
      jbinpartnerOfBin <- PERSIST.MX$hits[PERSIST.MX$hits[,"i"]==bin,"j"]
      ibinpartnerOfBin <- PERSIST.MX$hits[PERSIST.MX$hits[,"j"]==bin,"i"]
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
