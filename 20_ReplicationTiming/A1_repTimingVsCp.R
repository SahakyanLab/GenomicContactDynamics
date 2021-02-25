################################################################################
# Replication timing vs. Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/20_ReplicationTiming"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ReplicationTiming"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
id = "nontumor"
reptimePath = paste0(data.dir, "/replication_timing/out_clustering_combined/hg19/", 
                     id, "/RT_data_hg19.RData")
#reptimePath = "/Users/ltamon/Desktop/RT_data_hg19.RData"
out.dir = paste0(wk.dir, "/out_repTimingVsCp")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
HiC.res = 40000
Cp.v = 1:21
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
source(paste0(lib, "/makebp.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
rt <- readRDS(file=reptimePath)
if( all(chr.v%in%unique(rt$chroms)) ){
  rt <- rt[rt$chroms%in%chr.v,]
} else {
  stop("Not all chr.v in rt dataset.")
}

rt.res <- unique(diff( rt$starts[rt$chroms==chr.v[1]] ))

#-------------------Check rt dataset

# Filter out some bins
col.TF <- !colnames(rt)%in%c("chroms", "starts")
dropBin.TF <- rt$set.count < 3 | rt$point.count < 3 | rt$norm==1 | rt$norm.cA==1 | rt$norm.cB==1
rt[dropBin.TF,col.TF] <- NA
filteredOut.perc <- sum(dropBin.TF)/nrow(rt)*100
rm(col.TF, dropBin.TF)

# Check if rt data has all hg19 bins and in order 
chrLen.df <- read.delim(file=chrLenfile, header=TRUE)
rt$bin <- ceiling(rt$starts/rt.res)

for( chr in unique(rt$chroms) ){
  
  chr.len <- chrLen.df$length.bp[chrLen.df$chromosome==chr]
  tot.bin <- ceiling(chr.len/rt.res)
  
  chr.TF <- rt$chroms==chr

  if(
    !identical(as.numeric(rt$bin[chr.TF]), as.numeric(1:tot.bin))
  ){
    stop( paste0(chr, ": Checkpoint 1.") )
  }
  
  if( max(rt$starts[chr.TF]) > chrLen.df$length.bp[chrLen.df$chromosome==chr] ){
    stop( paste0(chr, ": Checkpoint 2.") )
  }
  
}

#-------------------Convert rt res to HiC res and obtain values per HiC res bin 
#(if res are not equal)

if( HiC.res%%rt.res==0 & rt.res<HiC.res ){
  
  rt$bin <- NULL
  rt$ind <- NA
  
  factr <- HiC.res/rt.res
  
  for( chr in unique(rt$chroms) ){
    
    chr.TF <- rt$chroms==chr
    ind <- rep(x=1:sum(chr.TF), each=factr)
    ind <- ind[1:sum(chr.TF)]
    rt[chr.TF, "ind"] <- ind
    rm(ind, chr.TF)
  }
  
  if( any(is.na(rt$ind)) ){
    stop("Checkpoint 4.")
  }
  
  # Calculate replication timing calculations per HiC res bin
  col.TF <- !colnames(rt)%in%c("chroms", "starts", "point.count", "ind")
  # NaN comes from taking mean of missing values
  # Wrong to add set.count because bins to be collapsed could share sets, so just take the mean
  rt1 <- aggregate(x=rt[,col.TF], by=list(chroms=rt$chroms, bin=rt$ind), FUN=mean, na.rm=TRUE)
  rt2 <- aggregate(x=rt$point.count, by=list(chroms=rt$chroms, bin=rt$ind), FUN=sum, na.rm=TRUE)
  
  if( !identical(rt1$chroms, rt2$chroms) ){
    stop("Checkpoint 5.")
  }
  if( !identical(rt1$bin, rt2$bin) ){
    stop("Checkpoint 6.")
  }
  
  # rt dataframe contains replication timing values per HiC res bin
  rt2$chroms <- NULL; rt2$bin <- NULL
  rt <- cbind(rt1, point.count=rt2$x)
  rm(rt1, rt2); gc()
  
} else if(HiC.res==rt.res){
  
  rt$starts <- rt$bin
  rt$bin <- NULL
  colnames(rt)[colnames(rt)=="starts"] <- "bin"
  print("Resolutions equal.", quote=FALSE)
  
} else {
  
  stop("Checkpoint 3.")

}
rownames(rt) <- paste0(rt$chroms, ".", rt$bin)

#-------------------Replication timing Vs. Cp
rt <- cbind.data.frame(
  rt, 
  matrix(data=NA, nrow=nrow(rt), ncol=length(Cp.v), dimnames=list(NULL, Cp.v))
)
if( any(duplicated(colnames(rt))) ){ stop("Checkpoint 7.") }

# Get Cp data of bin from BIN.MX
for(chr in chr.v){
  
  # BIN.MX
  load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
  BIN.MX <- BIN.MX[,!colnames(BIN.MX)%in%c("start","end")]
  rw.nme <- paste0(chr, ".", rownames(BIN.MX))
  
  # Check if order of bins in RT and BIN.MX is identical
  chr.TF <- rt$chroms==chr 
  if( !identical(rw.nme, rownames(rt)[chr.TF]) ){
    stop(paste0(chr, ": Checkpoint 8."))
  }
  rm(chr.TF)
  
  # Assign bin to Cp (involved in at least 1 contact for that Cp in any tissue)
  for( Cp in as.character(Cp.v) ){
    
    Cp.ind <- grep(x=colnames(BIN.MX), pattern=paste0("s_Cp_", Cp, "_ct"))
    if( length(Cp.ind)!=length(Cp.v) ){
      stop(paste0(chr, ": Checkpoint 9."))
    }
    rt[rw.nme, Cp] <- as.numeric(rowSums(x=BIN.MX[,Cp.ind])>0)
    rm(Cp.ind)
    
  }
  
  rm(BIN.MX, rw.nme); gc()
  print(paste0(chr, " done!"), quote=FALSE)
  
}
rownames(rt) <- rt$chroms <- rt$bin <- NULL

if( sum(is.na(rt[,as.character(Cp)]))>0 ){
  stop("Checkpoint 10.")
}

#col.TF <- !colnames(rt)%in%as.character(Cp.v)
#calc.v <- sort(colnames(rt)[col.TF])
calc.v <- c("mean", "mean.cA", "mean.cB", "median", "median.cA", "median.cB",   
            "norm", "norm.cA", "norm.cB", "sd", "sd.cA", "sd.cB", 
            "set.count", "point.count", "skew", "clust.p.value")
#rm(col.TF)
df <- reshape2::melt(data=rt, id.vars=calc.v)
# Choose only bins with Cp (only bins forming long-range contacts)
df <- df[df$value==1,]

#  Boxplot - replication timing value  vs. Cp
plot.id <- paste0(gcb, "_", id, "_rtres", rt.res, "_repTimingVsCp")
pdf(file=paste0(out.dir, "/", plot.id, ".pdf"), height=60, width=30)
par(mfrow=c(6,3))

for( calc in calc.v ){
  
  outline <- ifelse(calc=="set.count", TRUE, FALSE)
  makebp(df=df, x="variable", y=calc, xlab="Cp", ylab=calc, outline=outline,
         addjitter=FALSE, plot.id=paste0(plot.id, "_", calc, "_percFiltered", filteredOut.perc,
                                         "_totalbins", nrow(rt))
         )
  
}

dev.off()

# rm(list=ls()); gc()





