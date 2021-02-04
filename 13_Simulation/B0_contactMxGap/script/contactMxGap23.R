################################################################################
# Identify gaps in contact matrices so contacts can be masked appropriately
# and consistently across tissues.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
Cs.raw.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.disc.kmer.5.dir = CII.cont.kmer.5.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/polished/11_Constraints/out_group"
SIM.3.2.kmer.5.dir = paste0(wk.dir, "/sim_3.2")
SIM.4.2.kmer.5.dir = paste0(wk.dir, "/sim_4.2")
out.dir = paste0(wk.dir, "/out_contactMxGap")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chrX"
ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
              "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC"))
bin.len = 40000
# Metric name should match source directory name.
# 5 in "CII.disc.kmer.5" is the cutoff percentage for categorisation. disc means
# discrete (categorised CII), cont means continuouos (orig CII). 
# <CII/SIM>.<disc/cont>.<kmer/align>.<(0,100)>
metric = "Cs.norm"
out.id = "whole"
nCPU = 6L #~7G
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
library(GenomicRanges)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(wk.dir, "/lib/getContactDF.R"))
source(paste0(wk.dir, "/lib/filterContacts.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste(gcb, chr, out.id, metric, sep="_")
print(paste0(out.name, "..."), quote=FALSE)

eval(parse(text=paste0(
  'metric.dir <- ', metric, '.dir'
)))

genome <- read.table(file=chrLenfile, stringsAsFactors=FALSE, header=TRUE,
                     colClasses=c("character", "integer", "integer"))
mx.len <- ceiling(genome$length.bp[genome$chromosome==chr]/bin.len)
# Get reference to make sure df is upper triangle contacts only
if(mx.len!=genome$bins.40kb[genome$chromosome==chr]){ stop("mx.len wrong") } 
rm(genome)
temp <- expand.grid(1:mx.len, 1:mx.len)
temp <- temp[temp[,1]>temp[,2],]

toExport <- c("ct.v", "metric.dir", "metric", "gcb", "chr", "mx.len", "temp")
ct.v.len <- length(ct.v)
#### PARALLEL EXECUTION #########
X <- foreach(itr=isplitVector(1:ct.v.len, chunks=nCPU), 
             .inorder=FALSE, .combine="rbind.data.frame",
             .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
    
    ct <- ct.v[i]
    # Get upper triangle contacts
    df <- getContactDF(metric.dir=metric.dir, metric=metric, 
                       gcb=gcb, chr=chr, ct=ct, gap.range=NULL, 
                       incl.bin.x=NULL, incl.bin.y=NULL, 
                       mask.bin.x=NULL, mask.bin.y=NULL)
    if( nrow(df)!=((mx.len*mx.len)-mx.len)/2 ){
      stop("Wrong number of upper triangle contacts.") 
    }
    if( any(is.na(df$value)) & !grepl(x=metric, pattern="CII.|SIM.") ){
      stop("NA values in df.")
    }
    # Reorder
    df <- df[order(df$i, df$j),]
    # Converted to lower matrix df before comparing to reference
    colnames(df) <- c("j", "i", "value")
    # Make sure upper triangle contacts df contacts same orders as reference
    if( !identical(as.numeric(df$i), as.numeric(temp$Var1)) | 
        !identical(as.numeric(df$j), as.numeric(temp$Var2)) ){
      stop("Order of df wrong.")
    }
    
    # Add values to whole matrix
    MX <- matrix(data=NA, nrow=mx.len, ncol=mx.len)
    MX[ lower.tri(MX, diag=FALSE) ] <- as.numeric(df$value)
    MX <- t(MX)
    MX[ lower.tri(MX, diag=FALSE) ] <- as.numeric(df$value)
    rm(df); gc()
    diag(MX) <- 0
    
    row.sum <- rowSums(MX, na.rm=TRUE)
    row.ind <- which(row.sum==0) 
    col.sum <- colSums(MX, na.rm=TRUE)
    col.ind <- which(col.sum==0) 
    rm(MX, row.sum, col.sum); gc()
    
    y <- reduce(as(data.frame(chrom=chr, start=row.ind, end=row.ind), "GRanges"))
    y <- cbind.data.frame(as.data.frame(y), axis="y", stringsAsFactors=FALSE)
    x <- reduce(as(data.frame(chrom=chr, start=col.ind, end=col.ind), "GRanges"))
    x <- cbind.data.frame(as.data.frame(x), axis="x", stringsAsFactors=FALSE)
    # Combine
    x <- rbind.data.frame(x,y)
    rm(y, row.ind, col.ind)
    x <- cbind.data.frame(x, ct=ct, stringsAsFactors=FALSE)
    
    print(paste0(ct, " done!"), quote=FALSE)
    return(x)
    
  })
  return(do.call("rbind.data.frame", chunk))
  
}
### END OF PARALLEL EXECUTION ###
rm(temp, mx.len, metric.dir); gc()

X <- X[order(X$width, decreasing=TRUE),]
write.csv(X, file=paste0(out.dir, "/", out.name, "_contactMxGap.csv"), 
          row.names=TRUE)

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()

