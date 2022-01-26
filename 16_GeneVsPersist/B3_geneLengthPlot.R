################################################################################
# Title
# deva, R/3.5.0-newgcc, gcc/4.9.2
# deva, R/3.6.0-newgcc, gcc/4.9.2
# Mac, R/3.5.2, R/3.6.1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Set recommended global options

# Avoid left to right partial matching by $
#options(warnPartialMatchDollar=TRUE)

# Expands warnings
#options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/16_GeneVsPersist")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" 
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/3_AnnotationVsPersist")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")

data.dir = paste0(home.dir, "/Database")
len.dir = paste0(wk.dir, "/out_geneLength")
genelist.dir = paste0(wk.dir, "/z_ignore_git/out_anno_union")
out.dir = paste0(wk.dir, "/out_geneLengthPlot")

refseq = "ALL"
LTr.file = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19annoLTr_", refseq)
### OTHER SETTINGS #############################################################
anno.nme = "hg19anno"
gcb = "min2Mb"
nCPU = 1 # No. of cps
genelist.prefix = paste0(gcb, "_LTr")

#lentyp.v = c("len_full", "len_exon", "mean_exon", "num_exon",
#             "div_intronBYexon", "len_intron", "mean_intron", "num_intron")
#lentyp.v = c("len_full", "len_exon", "len_intron", 
#             "div_intronBYexon", "mean_exon", "mean_intron",
#             "div_intronBYexon", "num_exon", "num_intron")
#plot.id = "transcript_lengths"

lentyp.v = c("len_repeat_full", "len_repeat_exon", "len_repeat_intron",
             "fr_repeat_full", "fr_repeat_exon", "fr_repeat_intron")

plot.id = "repeat_values" 

# Plot parameters
repfree = F
outliers = F
dim.v = c(3,3)
plotOnly = F
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools) 
library(reshape2)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(wk.dir, "/lib/plotLENCPDF.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(file=paste0(len.dir, "/", anno.nme, "_ALL_annoLengths.RData"))
exclude.col <- c("name", "name2", "num_full", "num_full_repfree", "num_exon_repfree", 
                 "num_intron_repfree")
LEN.DF <- LEN.DF[,!colnames(LEN.DF)%in%exclude.col]

# Subset longest transcripts

LTr.df <- data.table::fread(file=LTr.file, header=T, data.table=F, stringsAsFactors=F)
LEN.DF <- LEN.DF[LEN.DF$uniqueID%in%LTr.df$uniqueID,]
rownames(LEN.DF) <- LEN.DF$uniqueID

rm(LTr.df, LTr.file, exclude.col)
gc()

# Get list of unique IDs per Cp

idPerCp <- readLines(con=paste0(genelist.dir, "/", genelist.prefix, "_", refseq,
                                "_uniqueID"))
header.ind <- grep(x=idPerCp, pattern=">all_genes_cp_", fixed=T)
cp.header <- idPerCp[header.ind]
idPerCp <- sapply(X=idPerCp[header.ind + 1], simplify=F, FUN=function(id.str){
  id <- strsplit(x=id.str, split=";", fixed=T)[[1]]
  return( as.character(id) )
})
cp.header <- gsub(x=cp.header, pattern=">all_genes_cp_|_end", replacement="")
names(idPerCp) <- cp.header

rm(cp.header, header.ind)

# Plot data and plot

lentyp.v.len <- length(lentyp.v)
out.id <- paste0(anno.nme, "_", refseq, "_", gcb)

pdf(file=paste0(out.dir, "/", out.id, "_", plot.id, ".pdf"),
    height=5*dim.v[[1]], width=5*dim.v[[2]])
par(mfrow=dim.v)

length.v <- colnames(LEN.DF)
for(lentyp in lentyp.v){
  
  len.v <- lentyp
  if( !grepl(x=lentyp, pattern="repeat") ){
    len.v <- c(lentyp, paste0(lentyp, "_repfree"))
    len.v <- intersect(length.v, len.v)
  }
  
  toExport <- c("idPerCp", "LEN.DF", "len.v")
  
  if(plotOnly==F){
    
    LENCP.DF <- foreach(itr=isplitVector(x=names(idPerCp), chunks=nCPU),
                        .inorder=T, .combine="rbind",
                        .export=toExport, .noexport=ls()[!ls()%in%toExport]
    ) %op% {
      
      chunk <- sapply(X=itr, simplify=FALSE, FUN=function(cp){
        
        mx <- cbind(Cp=as.numeric(cp), LEN.DF[idPerCp[[cp]], len.v])
        dimnames(mx)[[2]] <- c("cp", len.v)
        return(mx)
        
      })
      
      do.call("rbind", chunk)
      
    }
    
    rownames(LENCP.DF) <- NULL
    save(LENCP.DF, file=paste0(out.dir, "/", out.id, "_", lentyp, ".RData"))
    
  } else {
    load(file=paste0(out.dir, "/", out.id, "_", lentyp, ".RData"))
  }
  
  plotLENCPDF(DF=LENCP.DF, repfree=repfree, outliers=outliers)
  
  print(paste0(lentyp, " done!"), quote=F)
        
  rm(LENCP.DF, len.v, lentyp)
  gc()
  
}

dev.off()

# rm(list=ls()); gc()