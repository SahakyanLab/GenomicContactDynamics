################################################################################
# Make table containing (complete and repeat-free versions for each length):
# TRANSCRIPT.L <- transcript length
# N.EXONS <- exon count
# N.INTRONS <- intron count
# EXONS.L <- summed exon lengths
# INTRONS.L <- summed intron lengths
# MEAN.EXON.L <- mean single exon length 
# MEAN.INTRON.L <- mean single intron length 
# INTRONS.dev.EXONS <- total intron to exon ratio
# ANNOLENGTH.DF output contains data for all transcripts and data for longest 
# transcript can be subsetted using unique ID. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/16_GeneVsPersist"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = paste0(wk.dir, "/out_geneLength")
### OTHER SETTINGS #############################################################
# annotation file prefix
anno.nme = "hg19anno"
refseq = "ALL"
# For subtracting repeats
bedtoolsPATH="/Users/ltamon/prog/bedtools2/bin"
# Regenerate base ANNOLENGTH.DF (no repeat-free lengths)
regenerateBase=FALSE
# Regenerate beds subtracted with repeats?
regenerateBeds=FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)

diffNumString <- function(startString=annotable$exonStarts,
                          
                       endString=annotable$exonStarts,
                       splitChar=","){
  split.lst <- lapply(X=strsplit(x=c(startString, endString), split=splitChar),
                      FUN=as.numeric)
  return( split.lst[[2]]-split.lst[[1]] )
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(regenerateBase==TRUE){
  
  ANNOLENGTH.DF <- fread(file=paste0(annofile.dir, "/", anno.nme, "_", refseq), 
                         header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                         select=c("txStart", "txEnd", "exonStarts", "exonEnds",
                                  "exonCount", "name", "name2")) 
  
  ANNOLENGTH.DF <- within(ANNOLENGTH.DF, {
    
    # No +1 because UCSC anno tables are 0-based
    TRANSCRIPT.L <- txEnd-txStart
    txEnd <- NULL
    txStart <- NULL
    N.EXONS <- exonCount
    # Introns sandwiched by exons, due to the conventional process of splicing
    # tRNA processing is an exception
    N.INTRONS <- exonCount-1L
    exonCount <- NULL
    EXONS.L <- apply(X=ANNOLENGTH.DF[,c("exonStarts", "exonEnds")], MARGIN=1, 
                     FUN=function(row){
                       diff.v <- diffNumString(startString=row[1], endString=row[2], 
                                               splitChar=",") 
                       if(!all(diff.v>0)){
                         stop("Negative exon length.")
                       } else {
                         return(sum(diff.v, na.rm=TRUE))
                       }
                     })
    exonStarts <- NULL
    exonEnds <- NULL
    INTRONS.L <- TRANSCRIPT.L-EXONS.L
    MEAN.EXON.L <- EXONS.L/N.EXONS
    MEAN.INTRON.L <- INTRONS.L/N.INTRONS
    INTRONS.dev.EXONS <- INTRONS.L/EXONS.L
    
  })
  
  setnames(x=ANNOLENGTH.DF, old=c("name", "name2"), 
           new=c("RS.ACCESSION", "HUGO.SYMBOL"))
  
  # 4298th row, INTRONS.L == 0  
  if( !all(ANNOLENGTH.DF$TRANSCRIPT.L>0) | !all(ANNOLENGTH.DF$INTRONS.L>=0) ){
    stop("Negative values either in TRANSCRIPT.L or INTRONS.L")
  }
  
  # Reorder columns
  data.table::setcolorder(x=ANNOLENGTH.DF, neworder=c(1:2, 10:3))
  
  save(ANNOLENGTH.DF, file=paste0(out.dir, "/", anno.nme, "_", refseq, 
                                  "_annoLengths_base.RData"))
  
} else {
  load(file=paste0(out.dir, "/", anno.nme, "_", refseq, 
                   "_annoLengths_base.RData"))
}
### ADD REPEAT-FREE LENGTHS ####################################################

if(regenerateBeds==TRUE){
  
  repeat.bed <- fread(file=repmaskfile,
                      header=TRUE, data.table=FALSE, stringsAsFactors=FALSE, 
                      select=c("genoName", "genoStart", "genoEnd", 
                               "repName", "bin", "strand"),
                      col.names=c("chr", "start", "end", 
                                  "repName", "bin", "strand"))
  write.table(x=repeat.bed, file=paste0(out.dir, "/hg19_repeat.bed"), 
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  
  anno.bed <- fread(file=paste0(annofile.dir, "/", anno.nme, "_", refseq), 
                    header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                    select=c("chrom", "txStart", "txEnd", 
                             "name2", "uniqueID", "strand"),
                    col.names=c("chr", "start", "end", 
                                "name2", "uniqueID", "strand"))
  write.table(anno.bed, file=paste0(out.dir, "/", anno.nme, "_", refseq, "_tr.bed"), 
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  
  # Use subtractbed from bedtools to get Transcipt portions with no repeats
  # Inputs are the beds for each length
  # Output file = <genome.ver>anno_<refseq>_repFreeTr.bed
  # -s force same strandedness, -f min overlap 1bp (default)
  # Shell: subtractBed -s -f 1E-9 -a hg19anno_ALL_tr.bed -b hg19_repeat.bed > hg19anno_ALL_repFreeTr.bed
  # Shell: subtractBed -s -f 1E-9 -a hg19anno_ALL_exon.bed -b hg19_repeat.bed > hg19anno_ALL_repFreeEx.bed
  # Shell: subtractBed -s -f 1E-9 -a hg19anno_ALL_intron.bed -b hg19_repeat.bed > hg19anno_ALL_repFreeInt.bed
} 

# ANNOLENGTH.DF follows annotation file row order (only for refseq=ALL) 
# so uniqueID can be used as row number in this case 
fileid.v <- c("Tr", "Ex", "Int")
bednme.v <- c("TRANSCRIPT", "EXONS", "INTRONS")
v.len <- length(fileid.v)

for(i in 1:v.len){
  fileid <- fileid.v[i]
  bednme <- bednme.v[i]
  bed <- fread(file=paste0(out.dir, "/", anno.nme, "_", refseq, 
                           "_repFree", fileid, ".bed"), 
               header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
               col.names=c("chr", "start", "end", 
                           "name2", "uniqueID", "strand"))
  if(!all(bed$end>bed$start)){
    stop("Checkpoint 1.")
  }
  
  # Count number of exons and introns 
  #if(bednme!="TRANSCRIPT"){
  #  count.uniqueID <- stack( table(bed$uniqueID) )
  #  col1 <- paste0("N.", bednme, ".RF")
  #  ANNOLENGTH.DF[[col1]] <- NA
  #  ANNOLENGTH.DF[as.numeric( count.uniqueID$ind ), col1] <- count.uniqueID$values 
  #}
  
  bed$L.RF <- bed[,"end"]-bed[,"start"]
  agg <- aggregate(x=bed[,c("uniqueID", "L.RF")], by=list(bed$uniqueID), FUN="sum")[, c(1,3)]
  setnames(x=agg, old="Group.1", new="uniqueID")
  col1 <- paste0(bednme, ".L.RF")
  ANNOLENGTH.DF[,col1] <- NA
  ANNOLENGTH.DF[agg$uniqueID,col1] <- agg$L.RF
  
  col2 <- paste0("REP.PERC.", toupper(fileid), ".L")
  full <- paste0(bednme, ".L")
  ANNOLENGTH.DF[[col2]] <- 100*( ANNOLENGTH.DF[[full]]-ANNOLENGTH.DF[[paste0(full, ".RF")]] )/ANNOLENGTH.DF[[full]]
  
  rm(bed, agg, fileid, col1, col2, full); gc()
  
  print(bednme)
}

ANNOLENGTH.DF$INTRONS.dev.EXONS.RF <- ANNOLENGTH.DF$INTRONS.L.RF/ANNOLENGTH.DF$EXONS.L.RF

if( !all(unlist(ANNOLENGTH.DF[,-(1:2)])>=0, na.rm=TRUE) ){
  stop("Negative values in final ANNOLENGTH.DF.")
}

save(ANNOLENGTH.DF, file=paste0(out.dir, "/", anno.nme, "_", refseq, 
                                "_annoLengths.RData"))

# rm(list=ls())



  
