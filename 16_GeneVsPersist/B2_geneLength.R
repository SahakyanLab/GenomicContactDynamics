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
# LEN.DF output contains data for all transcripts and data for longest 
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
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/16_GeneVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = paste0(wk.dir, "/out_geneLength")
anno.dir = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/")
### OTHER SETTINGS #############################################################
# annotation file prefix
anno.nme = "hg19anno"
refseq = "ALL"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(file=paste0(out.dir, "/", anno.nme, "_ALL_annoLengths_base.RData"))

LEN.DF <- within(LEN.DF, {
  
  mean_exon <- len_exon/num_exon
  mean_exon_repfree <- len_exon_repfree/num_exon
  
  mean_intron <- len_intron/num_intron
  mean_intron_repfree <- len_intron_repfree/num_intron
  
  div_intronBYexon <- len_intron/len_exon
  div_intronBYexon_repfree <- len_intron_repfree/len_exon_repfree
  
})

LEN.DF <- LEN.DF[,sort(colnames(LEN.DF))]

anno.df <- data.table::fread(file=paste0(anno.dir, "/", anno.nme, "_", refseq),
                             data.table=F, header=T, stringsAsFactors=F)

if( identical(anno.df$uniqueID, LEN.DF$uniqueID) ){
  LEN.DF <- cbind.data.frame(name=anno.df$name, name2=anno.df$name2, LEN.DF)
} else {
  stop("anno.df and LEN.DF unique ids not identical.")
}

save(LEN.DF, file=paste0(out.dir, "/", anno.nme, "_ALL_annoLengths.RData"))

# rm(list=ls()); gc()