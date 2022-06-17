################################################################################
# Identify maximum Cp of contacts linking enhancer-target pairs
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/21_")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/21_")
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")

persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir  = paste0(wk.dir, "/out_scr")
anno.file = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno_ALL")
enhancer.file = paste0(data.dir, "/enhancers/ENdb_enhancer.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
bin.len = 40000
nCPU = 1
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
anno.df <- anno.df[,c("chrom", "txStart", "txEnd", "name2", "name", "uniqueID")]
anno.df$txStart <- anno.df$txStart + 1L

#
pair.df <- read.delim(file=enhancer.file, stringsAsFactors=F, header=T)
pair.df <- pair.df[pair.df$Reference_genome=="hg19",]
pair.df <- pair.df[,c("Chromosome", "Start_position", "End_position", 
                      "Enhancer_id", "Enhancer_symbol", "Target_gene")]

pair.df$Start_position <- pair.df$Start_position + 1L
pair.df$Target_gene <- enc2utf8(pair.df$Target_gene)
pair.df$Target_gene <- stringr::str_replace_all(string=pair.df$Target_gene,
                                                pattern=" ", replacement="") 

#
commaSepTarg <- strsplit(x=pair.df$Target_gene, split=",")
ind.v <- which(lengths(commaSepTarg) > 1)
tmp <- sapply(X=ind.v, simplify=F, FUN=function(ind){
  
  cbind(pair.df[ind, colnames(pair.df)!="Target_gene"], 
        `Target_gene`=commaSepTarg[[ind]], row.names=NULL)
  
})
tmp <- do.call("rbind", tmp)

pair.df <- pair.df[-ind.v,]
pair.df <- rbind(pair.df, tmp)

#
pair.df <- pair.df[toupper(pair.df$Target_gene) %in% toupper(anno.df$name2),]

#
pair.df <- merge(x=anno.df, y=pair.df, by.x="name2", by.y="Target_gene", all=F)
pair.df$Cp <- "None"

# a = gene, b = enhancer
p.coord.mx <- pair.df[,c("txStart", "txEnd", "Start_position", "End_position")]
p.coord.mx <- as.matrix(p.coord.mx)

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
  pair.df$Cp[is.chr][cp.mx[,1]] <- cp.mx[,2]
  
  rm(cp.mx, is.chr)
  
}

write.csv(pair.df, file=paste0(out.dir, "/", gcb, "_enhancer-target_maxCp.csv"), 
          row.names=F)

# rm(list=ls()); gc()