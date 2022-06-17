################################################################################
# Identify maximum Cp of contacts linking gene pairs from FunCoup networks
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Avoid left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/X1_")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/X1_")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")

persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir  = paste0(wk.dir, "/out_protein-protein1")
anno.file = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19annoLTr_ALL")
data.file = paste0(data.dir, "/FunCoup/FC5.0_H.sapiens_compact_GeneName")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
bin.len = 40000
nCPU = 4
out.name = "protein-protein"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/getConsensusCpOfPairs.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
anno.df <- read.delim(file=anno.file, stringsAsFactors=F, header=T)
anno.df <- anno.df[,c("chrom", "txStart", "txEnd", "name2")]
anno.df$txStart <- anno.df$txStart + 1L
anno.df$name2 <- toupper(anno.df$name2)
if( any(duplicated(anno.df$name2)) ){
  rm(pair.df, anno.df)
  stop("Duplicated gene names.")
}

#
pair.df <- read.delim(file=data.file, stringsAsFactors=F, header=T)
pair.df$X4.Gene1 <- toupper(pair.df$X4.Gene1)
pair.df$X5.Gene2 <- toupper(pair.df$X5.Gene2)

chr.anno <- setNames(object=anno.df$chrom, nm=anno.df$name2)
pair.df$chrom <- chr.anno[pair.df$X4.Gene1]
pair.df <- pair.df[( pair.df$X4.Gene1 %in% anno.df$name2 ) & 
                   ( pair.df$X5.Gene2 %in% anno.df$name2 ) & 
                     pair.df$chrom == unname(chr.anno[pair.df$X5.Gene2]), ]

txStarts <- setNames(object=anno.df$txStart, nm=anno.df$name2)
txEnds <- setNames(object=anno.df$txEnd, nm=anno.df$name2)
rm(anno.df)

# a = gene1, b = gene2
p.coord.mx <- cbind(txStarts[pair.df$X4.Gene1], txEnds[pair.df$X4.Gene1],
                    txStarts[pair.df$X5.Gene2], txEnds[pair.df$X5.Gene2])
dimnames(p.coord.mx)[[1]] <- NULL 

pair.df$Cp <- "None"

chr.v <- unique(pair.df$chrom)
chr.v <- intersect(paste0("chr", c(1:22, "X")), chr.v)
for(chr in chr.v){
  
  is.chr <- pair.df$chrom==chr
  
  # ij.mx
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  
  ij.mx <- cbind(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
  rm(PERSIST.MX)
  gc()
  
  #
  cp.mx <- getConsensusCpOfPairs(bin.len=bin.len, nCPU=nCPU, ij.mx=ij.mx,
                                 p.coord.mx=p.coord.mx[is.chr,], consensus.FUN=max)
  pair.df$Cp[is.chr][cp.mx[,"pair.ind"]] <- cp.mx[,"value"]
  
  rm(cp.mx, is.chr)
  
  print(paste0(chr, " done!"), quote=F)
  
}

write.csv(pair.df, file=paste0(out.dir, "/", gcb, "_", out.name, "_maxCp.csv"), row.names=F)

# rm(list=ls()); gc()