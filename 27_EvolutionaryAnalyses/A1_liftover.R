################################################################################
# Liftover Phylo-HMRF hg38 contacts to hg19
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) 
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
src.dir = paste0(data.dir, "/Phylo-HMRF/out_combine_splitPerChr") # 0-based
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/27_EvolutionaryAnalyses")
out.dir = paste0(wk.dir, "/out_liftover")
### OTHER SETTINGS #############################################################
src.id = "genome_state_Phylo-HMRF_mapping_contact50K_norm"
src.header = T
chrs = paste0("chr", c(17))
nCPU = 1
LOchain = "hg38ToHg19"
LOwidth.min.bp = 30000 
# Take only regions maintaining bin resolution (50000 in our case) after liftover 
# as well as output regions >= LOwidth.min.bp (latter is relevant for regions broken 
# down into multiple regions after conversion i.e. 1 long plus few short ones). 
# Ideally, choose value > half of original bin resolution to get 1 output region for 
# 1 input region.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)

library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))

library(rtracklayer)
source(paste0(lib, "/TrantoR_liftOver/GEN_liftOveR.R"))
source(paste0(lib, "/TrantoR_liftOver/GEN_liftOverLoadChain.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrs.len <- length(chrs)

#### PARALLEL EXECUTION #########

foreach(itr=isplitVector(1:chrs.len, chunks=nCPU), .inorder=F
        
) %op% {

  for(i in itr){
    
    chr <- chrs[[i]]
    out.id <- paste0(chr, "_", LOchain, "_LOwidth.min.bp", LOwidth.min.bp, "_", src.id)
    
    df <- fread(file=paste0(src.dir, "/", chr, "_", src.id, ".txt"), 
                header=src.header, stringsAsFactors=F, data.table=F)
    df[,2] <- df[,2] + 1L
    df[,4] <- df[,4] + 1L
    
    # Coordinate conversion
    
    orig.ij.num <- length(df[,1])
    iLO.df <- liftOveR(conversion=LOchain, 
                       space=as.character(df[,1]), 
                       start=as.numeric(df[,2]),
                       end=as.numeric(df[,3]),
                       strand=rep(x="*", times=orig.ij.num),
                       getchain=F,
                       rmchain=F,
                       returnGRangesList=F)
    iLO.df$seqnames <- as.character(iLO.df$seqnames)
    
    jLO.df <- liftOveR(conversion=LOchain, 
                       space=as.character(df[,1]), 
                       start=as.numeric(df[,4]),
                       end=as.numeric(df[,5]),
                       strand=rep(x="*", times=orig.ij.num),
                       getchain=F,
                       rmchain=F,
                       returnGRangesList=F)
    jLO.df$seqnames <- as.character(jLO.df$seqnames)
    
    pdf(file=paste0(out.dir, "/", out.id, ".pdf"), width=10, height=10)
    par(mfrow=c(2,1))
    
    plot(density(iLO.df$width), col="darkblue", cex.main=0.5, 
         main=paste0("iwidthvaluesBeforeLOwidth.min.bpfiltering\n", out.id))
    plot(density(jLO.df$width), col="darkred", cex.main=0.5,
         main=paste0("jwidthvaluesBeforeLOwidth.min.bpfiltering\n", out.id))
    
    dev.off()
    
    iLO.df <- iLO.df[iLO.df$width >= LOwidth.min.bp, ]
    jLO.df <- jLO.df[jLO.df$width >= LOwidth.min.bp, ]
    
    if( any( c(duplicated(iLO.df$group), duplicated(jLO.df$group)) ) ){
      rm(iLO.df, jLO.df, df)
      stop(paste0(chr, ": Multiple output regions for an input region. 
                  Consider increasing LOwidth.min.bp to > half of original resolution."))
    }
    
    # 
    
    df$group <- 1:orig.ij.num
    iLO.df <- merge(x=df[,"group", drop=F], y=iLO.df, by="group", all.x=T)
    jLO.df <- merge(x=df[,"group", drop=F], y=jLO.df, by="group", all.x=T)
    ijLO.df <- merge(x=iLO.df, y=jLO.df, by="group", all=T, suffixes=c(".i", ".j"))
    
    if( !identical(ijLO.df$group, df$group) ){
      rm(ijLO.df)
      stop(paste0(chr, ": Converted dataframe different order with original."))
    }
    
    # 
    
    i.len <- (ijLO.df$end.i - ijLO.df$start.i + 1L)
    j.len <- (ijLO.df$end.j - ijLO.df$start.j + 1L)
    
    if( any( min(c(i.len, j.len), na.rm=T) <= 0 ) ){
      rm(ijLO.df)
      stop(paste0(chr, ": Zero/Negative 1-based length/s of regions."))
    }
    
    # Populate with NAs contacts with overlapping regions
    
    range.len <- apply(ijLO.df[,c("start.i", "end.i", "start.j", "end.j")], 
                       MARGIN=1, 
                       FUN=function(x) diff(range(x, na.rm=F)) + 1L)
    
    is.ovrlapping <- range.len < (i.len + j.len) 
    is.ovrlapping[is.na(is.ovrlapping)] <- FALSE
    ijLO.df[is.ovrlapping,-1] <- NA 
    
    # Check chromosome correspondence
    
    is.validij <- !( is.na(ijLO.df$seqnames.i) | is.na(ijLO.df$seqnames.j) )
    if( !identical(ijLO.df$seqnames.i[is.validij], 
                   ijLO.df$seqnames.j[is.validij]) ){
      rm(ijLO.df)
      stop(paste0(chr, ": Inconsistent i and j chromosomes."))
    }
    
    # Switch non-overlapping, start i > end j
    
    is.idownstream <- (ijLO.df$start.i > ijLO.df$end.j)
    is.switch <- !is.ovrlapping & is.idownstream
    is.switch[is.na(is.switch)] <- FALSE
    cols.switch <- setdiff(colnames(iLO.df), "group")
    ijLO.df[is.switch,] <- ijLO.df[is.switch,
                                   c("group", 
                                     paste0(cols.switch, ".j"), 
                                     paste0(cols.switch, ".i"))]
    
    #
    ijLO.df$seqnames.j <- NULL
    ijLO.df$seqnames.i <- chr
    ijLO.df <- ijLO.df[, grepl(colnames(ijLO.df), pattern="seqnames.|start.|end.")]
    save(ijLO.df, file=paste0(out.dir, "/", out.id, ".RData"))
    
    print(paste0(chr, " done!"), quote=F)
    
  } 
  
}

### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()



