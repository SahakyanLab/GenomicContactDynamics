################################################################################
# Clean and reformat expression data for subsequent use.
# 1. Duplicated genes. There are very few duplicated genes in the expression data.
# I assumed them to correspond to multiple variants of the gene instead of being
# technical replicates because when comparing the values of some duplicated
# genes per tissue the values are so different e.g. 0.1 vs3 that it's dangerous
# to treat them as replicates and take the mean. We therefore just treated them
# as variant and took the sum of values per tissue. 
# 2. expr.cutoff argument. Expression below this cut-off is set to 0.
# 3. Histogram plots. Histogram plots were generated to show per gene, the number
# of tissues with data. Two histrograms were made, the first one uses unique genes
# in expression data and the other one uses only unique genes present in both 
# expression data and UCSC hg19 annotation table. 
# 4. chr. Using the UCSC annotation table, we added chr data to expression data.
# Genes with multiple transcripts in different chromosomes have their chromosomes
# written as a string separated by a period and genes not in the UCSC hg19 
# annotation table have NA or missing chr.
# 5. Final data. Final output is a csv file with chr, gene names and expression
# value for each tissue with expr.cutoff applied. The final table contains
# all genes in the original expression except that duplicated genes were resolved
# as described above. 
# 6. Expression datasets. I worked with two sets of baseline expression data 
# from EMBL-EBI expression atlas. data1 has 27 tissues while data2 has 53 tissues. 
# Values are in TPM. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir = "/Users/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
out.dir = paste0(wk.dir, "/out_cleanExprData")
gene.id = "LTr_ALL" #"LTr_ALL" "_ALL"
annofilePath = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
src.id = "data1"
expr.cutoff = 0 # Could use for both TPM and FPKM, expression value less than
# cut-off converted to 0; NA in table means that there is no data available
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0("expr_", src.id, "_cutoff", expr.cutoff, "_", gene.id)

# Path to protein-coding gene expression data (bulk RNA-seq)
if(src.id=="data1"){
  
  # 43706 unique gene names
  exprDataPath = paste0(data.dir, "/expression/RNA-seq_mRNA/data1/E-MTAB-1733-query-results.tpms.tsv")
  
} else if(src.id=="data2"){
  
  # 46704 unique gene names
  exprDataPath = paste0(data.dir, "/expression/RNA-seq_mRNA/data2/E-MTAB-5214-query-results.tpms.tsv")
  
} else {
  stop("exprDataPath for given src.id argument not supplied.")
}

expr.df <- read.delim(file=exprDataPath, skip=4, header=T, stringsAsFactors=F)

if( any(duplicated(expr.df$Gene.ID)) ){
  stop("Duplicated gene ids.")
}

expr.df <- expr.df[,colnames(expr.df)!="Gene.ID"]

#-------------------

# There duplicated gene names and we assued that they correspond to transcript
# variants so we took the sum of values per tissue. 
# 17 (data 1) 7 (data 2) duplicated genes
dup.genes <- unique(expr.df[duplicated(expr.df$Gene.Name), "Gene.Name"])
tisscol.TF <- colnames(expr.df)!="Gene.Name"

temp <- sapply(X=dup.genes, simplify=F, FUN=function(dup.gene){
  
  x <- apply(X=expr.df[expr.df$Gene.Name==dup.gene, tisscol.TF], 
             MARGIN=2, FUN=sum, na.rm=T)
  x[is.nan(x)] <- NA
  return(x)
  
})
temp <- do.call("rbind", temp)

expr.df <- expr.df[!duplicated(expr.df$Gene.Name),]
rownames(expr.df) <- expr.df$Gene.Name
expr.df[dup.genes, tisscol.TF] <- temp
rm(temp, dup.genes)
gc()
if( any(duplicated(expr.df$Gene.Name)) ){
  stop("Duplicated genes in filtered expr.df.")
}

#-------------------

# Apply expression value cut-off
for(i in 2:ncol(expr.df)){
  expr.df[expr.df[,i]<expr.cutoff & !is.na(expr.df[,i]) ,i] <- 0
}

#-------------------

# Check if genes present have data for at least one tissue
NtissuesWD <- rowSums(!is.na(expr.df[,tisscol.TF]))
if( any(NtissuesWD==0) ){
  stop("There are gene/s with no expression data in all tissues.")
}

pdf(file=paste0(out.dir, "/", out.name, "_datasetinfo.pdf"),
    height=10, width=20)
par(mfrow=c(1,2))

hist(NtissuesWD, cex.main=0.5, 
     main=paste0(out.name, 
                 "_TotalGenesInHist=NumGenesInNoDupsExprData", 
                 length(expr.df$Gene.Name), 
                 "_Ntissues", sum(tisscol.TF), 
                 "_minNtissuesWD", min(NtissuesWD), 
                 "_maxNtissuesWD", max(NtissuesWD), 
                 "\nEachValueInHistIsNumberOfTissuesWithDataPerGene"))

#-------------------

# Determine chr location of genes
anno.df <- read.delim(file=annofilePath, stringsAsFactors=F, header=T)[,c("chrom", "name2")]
common.genes <- unique(intersect(anno.df$name2, expr.df$Gene.Name))
NtissuesWD <- NtissuesWD[common.genes]

hist(NtissuesWD, cex.main=0.5, 
     main=paste0(out.name, 
                 "_TotalGenesInHist=NumUniqueCommonGenesOfUCSCAndNoDupsExprData", 
                 # intersect() discards duplicated values but used unique
                 # for peace of mind :)
                 length(NtissuesWD),
                 "_NuniqueUCSCGenes",
                 length(unique(anno.df$name2)),
                 "\nNtissues", sum(tisscol.TF),
                 "_minNtissuesWD", min(NtissuesWD), 
                 "_maxNtissuesWD", max(NtissuesWD), 
                 "_EachValueInHistIsNumberOfTissuesWithDataPerGene"))

dev.off()

#-------------------

# Add chr data to expr.df. For genes with multiple transcripts located 
# in different chr (not the case when using longest transcript per gene)
# store chr as period separated string. Genes in expr.df not in anno.df
# will have NA or missing chr.
expr.df <- cbind(chr=NA, expr.df)
expr.df$chr <- sapply(X=expr.df$Gene.Name, simplify=T, FUN=function(expr.gene){
  
  chr.v <- sort(unique(anno.df[anno.df$name2==expr.gene, "chrom"]))
  if(length(chr.v)>0){
    return( paste(paste0(chr.v, "."), collapse="") )
  } else {
    return(NA)
  }
  
})
rm(anno.df)
gc()

#-------------------

expr.cutoff <- gsub(x=expr.cutoff, pattern=".", replacement="", fixed=T)
write.csv(expr.df, file=paste0(out.dir, "/", out.name, ".csv"),
          quote=F, row.names=F)

# rm(list=ls()); gc()



