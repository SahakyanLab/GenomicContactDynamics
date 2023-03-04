################################################################################
# Per element, calculate stats pertaining to genome coverage and range of site
# length. Append also copy number calculated previously.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
repmask.file = paste0(data.dir, "/ucsc_tables/hsa_RepeatMasker/RepeatMasker_hg19")
copyNum.file = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/17_RepeatAge/out_repeatContent/hg19repeats_copyNumber.RData"
out.dir = paste0(wk.dir, "/out_site_coverage_length")
### OTHER SETTINGS #############################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(GenomicRanges)
library(Rmisc)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
repmasker.df <- fread(file=repmask.file, header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)
elements <- unique(repmasker.df$repName) # REMOVE
elm.len <- length(elements)

#all(repmasker.df$genoEnd > repmasker.df$genoStart)

mx <- matrix(data=NA, nrow=elm.len, ncol=6)

for(i in 1:elm.len){
  
  elm <- elements[[i]]
  elm.df <- repmasker.df[repmasker.df$repName==elm,]
  elm.gr <- makeGRangesFromDataFrame(df=elm.df, keep.extra.columns=F,
                                     ignore.strand=T, seqnames.field="genoName",
                                     start.field="genoStart", end.field="genoEnd",
                                     strand.field=NULL, starts.in.df.are.0based=T)
  
  # Genome coverage
  
  elm.lens <- width(elm.gr)
  mx[i,1] <- sum(elm.lens)
  mx[i,2] <- sum(width(GenomicRanges::reduce(elm.gr))) 
  
  # Length stats
  mx[i,3:6] <- c(min(elm.lens), max(elm.lens), mean(elm.lens), median(elm.lens))
  
  rm(elm, elm.df, elm.gr, elm.lens)
  
  message(paste0(i, ": done!"))
  
  #return(
  #  as.list(c(total.bp=total.bp, reduced.bp=reduced.bp,
  #            min.len.bp=min(elm.lens),
  #            max.len.bp=max(elm.lens),
  #            mean.len.bp=mean(elm.lens),
  #            med.len.bp=median(elm.lens)))
  #)
  
}

dimnames(mx) <- list(elements, c("total.bp", "reduced.bp", 
                                 "min.len.bp", "max.len.bp",
                                 "mean.len.bp", "median.len.bp"))

# Add copynumber data

load(copyNum.file)
rownames(COPYNUM$repName) <- COPYNUM$repName$repName
mx <- cbind(mx, copynum=COPYNUM$repName[ elements, "copyNumber"]) # REMOVE
mx <- mx[sort(dimnames(mx)[[1]], decreasing=F), ]

save(mx, file=paste0(out.dir, "/hg19_site_coverage_length.RData"))

# rm(list=ls()); gc()

