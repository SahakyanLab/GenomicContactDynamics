################################################################################
# Make of plot of transcript's length and count of contacts overlapping with it
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    annofile.dir = "/t1-data/user/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
trcount.dir = paste0(wk.dir, "/out_HiCTrExploration")
out.dir = paste0(wk.dir, "/out_Transcript_LenVsContOlap")
### OTHER SETTINGS #############################################################
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
gcb.v = "min2Mb" #min05Mb
chr = "chrALL"
refseq = #c("ALL", "NM", "NR")
anno.nme = "hg19anno"
# Max is total number of contacts
nCPU = 2L
# TRCOUNT.MX suffix
data.suffix = "ContactCountPerTr_EqTrLen"
out.name = "TRMEANCP_EqTrLen"
plotOnly = FALSE
plotOnlyLtr = FALSE
  datapoints = 1000L
### FUNCTION ###################################################################
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools) #isplitVector
library(ggplot2)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
annotable <- fread(file=paste0(annofile.dir, "/", anno.nme, "_", refseq), 
                    header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                    select=c("uniqueID", "txStart", "txEnd"))
uniqueID.v <- annotable$uniqueID

for(gcb in gcb.v){
  
  if(plotOnly==FALSE){
    
    suffix <- paste0(gcb, "_", refseq)
    
    # Load TRCOUNT.MX
    load(file=paste0(trcount.dir, "/", chr, "_", data.suffix, "_", suffix, 
                     ".RData"))
    
    #### FOREACH EXECUTION #########
    
    TRCPSUM <- foreach(itr=isplitVector(1:21, chunks=nCPU), 
                              .inorder=TRUE, .combine="+",
                              .export="TRCOUNT.MX", 
                          .noexport=ls()[!ls()%in%"TRCOUNT.MX"]
    ) %op% {
      
      chunk <- sapply(X=itr, simplify=FALSE, FUN=function(ntis){
        print(ntis)
        TRCOUNT.MX[,as.character(ntis)]*ntis
      })
      
      rowSums(do.call("cbind", chunk))
    }
    
    ### END OF FOREACH EXECUTION ###
    
    TRCOUNTSUM <- rowSums(TRCOUNT.MX)
    TRMEANCP.MX <- cbind(uniqueID=as.numeric(rownames(TRCOUNT.MX)),
                         meanCp=TRCPSUM/TRCOUNTSUM,
                         numContOlap=TRCOUNTSUM)
    rownames(TRMEANCP.MX) <- NULL
    
    rm(TRCPSUM, TRCOUNTSUM); gc()
    save(TRMEANCP.MX, file=paste0(out.dir, "/", chr, "_", suffix, 
                            "_", out.name, ".RData"))
    
  } else {
    
    load(file=paste0(out.dir, "/", chr, "_", suffix, "_", out.name, 
                     ".RData"))
  }
  
  if(plotOnlyLtr==TRUE){
    # PLOT (only for longest transcript)
    annotable <- fread(file=paste0(annofile.dir, "/", anno.nme, "LTr_", refseq), 
                       header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                       select=c("uniqueID", "txStart", "txEnd"))
    annotable$txSize <- annotable$txEnd - annotable$txStart
    uniqueID.v <- as.numeric(annotable$uniqueID)
    # sample() has a bug, if length(x)==1, it will use the vector 1:size
    # for the sampling
    uniqueID.v <- sample(x=uniqueID.v, size=datapoints)
    rm(annotable); gc()
    
    df <- data.frame(numContOlap=TRMEANCP.MX[TRMEANCP.MX[,"uniqueID"]%in%uniqueID.v,
                                                 "numContOlap"],
                     txSize=annotable[annotable$uniqueID%in%uniqueID.v, "txSize"], 
                     stringsAsFactors=FALSE)
  } else {
    # All transcripts
    annotable <- fread(file=paste0(annofile.dir, "/", anno.nme, "_ALL"), 
                       header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                       select=c("uniqueID", "txStart", "txEnd"))
    annotable$txSize <- annotable$txEnd - annotable$txStart
    df <- data.frame(numContOlap=TRMEANCP.MX[,"numContOlap"],
                     txSize=annotable[TRMEANCP.MX[,"uniqueID"], "txSize"], 
                     stringsAsFactors=FALSE)
    rm(annotable)
  }
    
  ggplot(data=df, aes(x=numContOlap, y=txSize)) +
    geom_point(color="gray59", shape=1) +
    stat_smooth(geom="smooth", method="lm", formula=y~x, se=TRUE, 
                n=80, span=0.75, level=0.95, colour="darkred", size=3) +
    ggtitle(label=paste0(chr, "_", suffix, "_Tr_LenVsContOlap_N=", nrow(df))) +
    ylab(label=expression("L"^"tr")) +
    xlab(label="No. of contacts") +
    bgr2 
  
  ggsave(filename=paste0(out.dir, "/", chr, "_", suffix, 
                         "_Tr_LenVsContOlap.pdf"), units="in",
         width=13, height=12)
  
} # gcb.v for loop end

# rm(list=ls())


  



