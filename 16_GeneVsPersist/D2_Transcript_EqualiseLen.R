################################################################################
# Equalise length of transcripts and make an annotation table of 
# chr-unique ID-HUGO symbol-newTxStart-newTxEnd (saved in Database)
# Equalise based on a reference length (Lref):
# a. mean of trancript lengths (ave)
# b. 2 SDs from mean of transcript lengths (m2sd)
# Also make a plot of transcript lengths with the Lref indicated 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    data.dir = "/Users/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
annofile.dir = paste0(data.dir, "/ucsc_tables/hsa_geneAnno")
out.dir = paste0(wk.dir, "/out_Transcript_EqualiseLen")
### OTHER SETTINGS #############################################################
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
refseq = "ALL" 
anno.nme = "hg19anno"
Lref.meth = "average" #c("average", "m2sd")
out.name = "EqTrLen_ave"
#out.name = "EqTrLen_m2sd"
trDistrPlot = TRUE
### FUNCTION ###################################################################
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# For calculation of Lref, consider only the longest transcript per gene
# because there's bias for longer genes to have more isoforms
LTrAnnot <- fread(file=paste0(annofile.dir, "/", anno.nme, "LTr_", refseq), 
                  header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                  select=c("txStart", "txEnd"))

LTrAnnot$txSize <- LTrAnnot$txEnd-LTrAnnot$txStart

if(Lref.meth=="average"){
  Lref <- ceiling( mean(LTrAnnot$txSize, na.rm=TRUE) )
  print("Lref is average gene length.", quote=FALSE)
} else if(Lref.meth=="m2sd"){
  mean.val <- mean(LTrAnnot$txSize, na.rm=TRUE)
  sdev.val <- sd(LTrAnnot$txSize, na.rm=TRUE)
  Lref <- as.integer( ceiling(mean.val+(2*sdev.val)) ) 
  print("Lref is two standard deviations from mean gene length.", quote=FALSE)
} else {
  print("Invalid input for Lref.method.")
}
rm(LTrAnnot, mean.val, sdev.val)

annotable <- fread(file=paste0(annofile.dir, "/", anno.nme, "_", refseq), 
                   header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                   select=c("chrom", "name2", "uniqueID", "txStart", "txEnd"))
annotable$txSize <- annotable$txEnd-annotable$txStart

if(trDistrPlot==TRUE){
  
  ggplot(data=as.data.frame(annotable$txSize), 
         aes(x=annotable$txSize, y=..count..)) +
    geom_density(stat="density", position="identity", na.rm=FALSE, 
                 bw="nrd0", adjust=1, kernel="gaussian", 
                 n=512, trim=FALSE, fill="deepskyblue3") +
    geom_vline( linetype="dashed", colour="black", size=0.7, 
                aes(xintercept=Lref) ) +
    #annotate("text", x=Lref+5e5, y=0.8, size=10,
    #         label = "paste(L ^ ref, \" = 280165\")", parse=TRUE) +
    ylab(label="Density (count)") +
    scale_x_continuous(name=expression("L"^"tr")) +
    ggtitle(paste0(anno.nme, out.name, "_TrDist_", refseq, "_Lref=", Lref)) + 
    bgr2 
  #breaks=c(0, Lref, 1e6, 2e6)
  ggsave(filename=paste0(out.dir, "/", anno.nme, out.name,
                         "_TrDist_", refseq, ".pdf"),
         units="in", width=9, height=9)
  
} # trDistrPlot

# Equalize length of transcripts to Lref


#annotable$txStart <- annotable$txStart+len.ext
#annotable$txEnd <- annotable$txEnd+len.ext
#annotable$txSize <- annotable$txEnd-annotable$txStart

# Chromosome lengths
chr.len <- fread(file=paste0(data.dir, "/Hsa_GRCh37_73_chr_info.txt"), 
                 header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                 select=c("chromosome", "length.bp"), 
                 col.names=c("chrom", "length"))
chr.v <-  unique(annotable$chrom)

#### FOREACH EXECUTION #########

annotable.eq <- foreach(chr=chr.v, .combine="rbind", .inorder=TRUE
        
) %do% {
  
  annot.sub <- annotable[annotable$chrom==chr,]
  chrlen <- chr.len[chr.len$chrom==chr, "length"]
  len.ext <- (Lref-annot.sub$txSize)/2
  
  mx <- cbind(
    uniqueID=annot.sub$uniqueID,
    ext=len.ext,
    txStart=annot.sub$txStart,
    txEnd=annot.sub$txEnd,
    start=as.integer( trunc(annot.sub$txStart-len.ext) ),  
    end=as.integer( trunc(annot.sub$txEnd+len.ext ) )
  )
  
  # Impossible to have ranges falling both in Cases 1 and 2
  
  # Case 1 - negative start points after extension
  ind1 <- which(mx[,"start"]<0)
  if( length(ind1)!=0 ){
    
    abs.val <- abs(mx[ind1, "start"])
    adj.val <- max(abs.val)-abs.val
    mx[ind1,"start"] <- 0 + adj.val
    mx[ind1,"end"] <- Lref + adj.val
    if( any(mx[ind1,"start"]<0) ){
      stop("Checkpoint 1")
    }
  }
  
  # Case 2 - end points exceeding chr length after extension
  
  ind2 <- which(mx[,"end"]>chrlen)
  if( length(ind2)!=0 ){
    adj.val <- max(mx[ind2, "end"], na.rm=TRUE)-mx[ind2,"end"]
    mx[ind2,"end"] <- chrlen-adj.val
    mx[ind2,"start"] <- (chrlen-Lref)-adj.val
    if( any(mx[ind2,"end"]>chrlen) ){
      stop("Checkpoint 2")
    }
  }
  
  print( paste0( chr, ":", length(unique(c(ind1, ind2))) ) )
  mx <- cbind(mx, size=trunc( mx[,"end"]-mx[,"start"] ) )
  
  # Check if adjusted ranges contain the transcript for extension>0
  ind.negext <- which(mx[,"ext"]<0)
  if(any( mx[-ind.negext,"txStart"]<mx[-ind.negext,"start"] & 
          mx[-ind.negext,"txEnd"]>mx[-ind.negext,"end"] &
          unique(mx[,"size"]!=Lref) ) 
  ){
    stop("Checkpoint 3.")
  }
  
  # Check if adjusted ranges within the transcript for extension<0
  if(any( mx[ind.negext,"txStart"]>mx[ind.negext,"start"] & 
          mx[ind.negext,"txEnd"]<mx[ind.negext,"end"] ) 
  ){
    stop("Checkpoint 4.")
  }
  
  mx <- mx[,-c(2:4,7)]
  dimnames(mx)[[2]] <- list("uniqueID", "txStart", "txEnd")
  return(mx)
}
  
annotable.eq <- merge(y=annotable.eq, x=annotable[, c("chrom", "uniqueID", "name2")],
                      by="uniqueID", all.x=TRUE)
if(!identical(annotable.eq$uniqueID, annotable$uniqueID)){
  stop("Checkpoint 5.")
}
#ordr <- match(x=annotable.eq$uniqueID, table=annotable$uniqueID)
#annotable.eq <- annotable.eq[order(ordr),]

### END OF FOREACH EXECUTION ###

write.table(x=annotable.eq, file=paste0(annofile.dir, "/", anno.nme, out.name, 
                                        "_", refseq),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# rm(list=ls())

