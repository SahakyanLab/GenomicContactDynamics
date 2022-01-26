################################################################################
# Identify longest transcript per gene. In case of ties, pick first entry - order
# based on original UCSC table but preferring coding over non-coding transcripts. 
# (962/24910 ~ 3.9% instances that we have to pick NM over NR)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/16_GeneVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
out.dir = paste0(wk.dir, "/out_LTr")
### OTHER SETTINGS #############################################################
refseq.v = c("ALL", "NM", "NR")
anno.nme = "hg19anno"
### FUNCTION ###################################################################
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(refseq in refseq.v){
  
  annotable <- fread(file=paste0(annofile.dir, "/", anno.nme, "_", refseq), 
                     header=T, data.table=F, stringsAsFactors=F)
  
  # For selecting coding/NM LTr over non-coding/NR, if applicable
  annotable$name <- gsub(x=annotable$name, replacement="", pattern="[_].*$")
  if( any(!unique(annotable$name)%in%c("NM", "NR")) ){
    stop(paste0(refseq, ": Error in renaming Refseq ID."))
  }
  
  if( !is.numeric(annotable$uniqueID) ){
    stop(paste0(refseq, ": Unique ID not numeric."))
  }
  
  # 0-based coordinates
  annotable$txSize <- annotable$txEnd-annotable$txStart
  if( any(annotable$txSize<=0) ){
    stop(paste0(refseq, ": Invalid transcript length."))
  }
  
  # Verify that NR ids non-coding
  annotable$cdsSize <- annotable$cdsEnd-annotable$cdsStart
  if( any(annotable$cdsSize[annotable$name=="NR"]!=0) | any(annotable$cdsSize[annotable$name=="NM"]<=0) ){
    stop(paste0(refseq, ": Error in NM and NR classification."))
  }
 
  LtrPerHugo <- by(data=annotable[,c("name", "uniqueID", "txSize", "cdsSize", "chrom", "txStart", "txEnd")],
                   INDICES=annotable$name2, simplify=F, # output is a list 
                   FUN=function(df){
                     
                     df <- df[df$txSize==max(df$txSize, na.rm=F),]
                     
                     df$numLTr <- 1
                     df$numchromLTr <- 1
                     df$multcoordLTr.TF <- 0
                     df$NMNRLTr.TF <- 0
                     
                     # Count problematic genes
                     
                     LTr.len <- length(df$chrom)
                     if(LTr.len>1){
                       
                       df$numLTr <- LTr.len
                       
                       uniqchrom.len <- length(unique(df$chrom))
                       if(uniqchrom.len>1){
                         df$numchromLTr <- uniqchrom.len
                       }
                       
                       if( length(unique(df$txStart))>1 | length(unique(df$txEnd))>1 ){
                         df$multcoordLTr.TF <- 1
                       }
                       
                       if( length(unique(df$name))>1 ){
                         
                         df$NMNRLTr.TF <- 1
                         df <- df[df$name=="NM",]
                         
                       } 
                       
                     }
                     
                     # Pick first transcript from top equal-length ones (as in 
                     # Sahakyan et al. 2016 BMC Genomics) so the one that will be picked
                     # depends on the order of the UCSC table but preferring coding over
                     # non-coding ones.
                     return( df[1,c("uniqueID", "numLTr", "numchromLTr", "multcoordLTr.TF", "NMNRLTr.TF")] )
                     
                   }) 
  
  LtrPerHugo <- do.call("rbind", LtrPerHugo)
  
  LtrPerHugo$multchromLTr.TF <- as.numeric(LtrPerHugo$numchromLTr>1)
  LtrPerHugo$multLTr.TF <- as.numeric(LtrPerHugo$numLTr>1)
  
  #annotable <- annotable[annotable$uniqueID%in%LtrPerHugo,]
  annotable <- annotable[annotable$uniqueID%in%as.numeric(as.character(LtrPerHugo$uniqueID)),]
  #write.table(x=annotable, file=paste0(annofile.dir, "/", anno.nme, "LTr_", refseq),
  #            sep="\t", quote=F, col.names=T, row.names=F)
  
  # REMOVE
  LTr.current <- fread(file=paste0(annofile.dir, "/", anno.nme, "LTr_", refseq), 
                       header=T, data.table=F, stringsAsFactors=F)
  id.curr <- LTr.current$uniqueID
  if( !identical(id.curr, annotable$uniqueID) ){
    stop(paste0(refseq, ": New code output different from current LTr files."))
  }
  
  if( any(duplicated(annotable$uniqueID)) ){
    stop(paste0(refseq, ": Duplicated LTr."))
  }
  
  LTRPERHUGO <- list(Ltr=LtrPerHugo,
                     problematic=table(LtrPerHugo[,c("multchromLTr.TF", "multcoordLTr.TF", "NMNRLTr.TF", "multLTr.TF")])
                     )
  if( sum(as.numeric(LTRPERHUGO$problematic))!=length(LTRPERHUGO$Ltr$uniqueID) ){
    stop(paste0(refseq, ": Checkpoint 1."))
  } else {
    save(LTRPERHUGO, file=paste0(out.dir, "/", anno.nme, "LTr_", refseq, "_data.RData"))
  }
  
  print(paste0(refseq, " done!"), quote=F)
  
  rm(annotable, LtrPerHugo, LTRPERHUGO, refseq)
  gc()
  
}

# rm(list=ls()); gc()


