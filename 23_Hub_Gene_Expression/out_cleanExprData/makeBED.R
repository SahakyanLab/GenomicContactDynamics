################################################################################
# BED file of cross-tissue gene expression value. From annotation file containing
# all transcripts, take those with expression data  only. For each entry assign the
# cross-tissue mean or median expression value. In the final BED, transcripts are
# individually listed to match shiny app gene information but since the expression 
# data is not differentiated by transcript, all transcripts of a gene will have 
# same expression value. 
################################################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/23_Hub_Gene_Expression/out_cleanExprData"
anno.file = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno/hg19anno_ALL"
out.dir = paste0(wk.dir, "/out_makeBED")

data.id = "data2" # data1 | data2
consensus.funx = "median" # mean | median
src.file = paste0("expr_", data.id, "_cutoff0_LTr_ALL")
out.id = paste0("hg19anno_ALL", src.file, "_crosstissue", consensus.funx)

################################################################################
anno.df <- read.delim(file=anno.file, header=T, stringsAsFactors=F)[,c("chrom", "txStart", "txEnd", "name2", "name")]
anno.df$NAME2 <- toupper(anno.df$name2)

df <- read.csv(file=paste0(wk.dir, "/", src.file, ".csv"), header=T, stringsAsFactors=F)
df <- cbind.data.frame(GENE.NAME=toupper(df$Gene.Name), df)
df <- df[ df$GENE.NAME %in% anno.df$NAME2, ]
anno.df <- anno.df[ anno.df$NAME2 %in% df$GENE.NAME, ]

if( any( lengths( strsplit(x=df$chr, split=".", fixed=T) ) != 1 ) ){
  stop("Multiple chr for gene/s.")
  rm(df)
}

if( any(duplicated(df$Gene.Name)) ){
  stop("Duplicated gene names.")
  rm(df)
}

#
ntiss.val <- apply( X=df[,-(1:3)], MARGIN=1, FUN=function(rw) sum(!is.na(rw)) )
eval(parse(text=paste0(
  'consensus.val <- apply( X=df[,-(1:3)], MARGIN=1, FUN=', consensus.funx,', na.rm=T)'
)))

#
df <- cbind.data.frame(GENE.NAME=df$GENE.NAME, 
                       ntiss.val_consensus.val=paste0(ntiss.val, "_", consensus.val), 
                       consensus.val=consensus.val)
bed <- merge(x=anno.df, y=df, by.x="NAME2", by.y="GENE.NAME", all=T)
bed <- bed[,c("chrom", "txStart", "txEnd", "ntiss.val_consensus.val", "consensus.val", "name2", "name")]

write.table(bed, file=paste0(out.dir, "/", out.id, ".bed"), sep="\t", row.names=F, col.names=F, quote=F)
################################################################################

# rm(list=ls()); gc()