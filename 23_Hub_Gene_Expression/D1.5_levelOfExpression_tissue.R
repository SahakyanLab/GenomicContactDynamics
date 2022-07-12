################################################################################
# Plot proportion of genes with no data, no, low, medium and high expression
# per tissue per expression dataset. Proportion is relative to total unique 
# genes in the UCSC hg19 annotation table.
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
exprData.dir = paste0(wk.dir, "/out_cleanExprData")
out.dir = paste0(wk.dir, "/out_levelOfExpression_tissue")

gene.id = "LTr_ALL" #"ALL"
annofilePath = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/hg19anno", gene.id)
### OTHER SETTINGS #############################################################
src.id = "data1"
expr.cutoff = 0 # Could use for both TPM and FPKM, expression value less than
# cut-off converted to 0
col.v = rev(c(ND="gray50", NE="gray70", LE="#2171B5", ME="#9ECAE1", HE="#A50F15"))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0("expr_", src.id ,"_cutoff0_", gene.id)

exprDataPath <- paste0(exprData.dir, "/expr_", src.id ,"_cutoff0_", gene.id, ".csv")
# Rows with NAs in all tissues were already removed
expr.df <- read.csv(file=exprDataPath, header=T, stringsAsFactors=F)
expr.df$chr <- NULL

anno.df <- read.delim(file=annofilePath, stringsAsFactors=F, header=T)[,"name2", drop=F]
expr.df <- expr.df[expr.df$Gene.Name%in%anno.df$name2,]

# Append to expr.df missing genes from anno.df

toaddgenes <- setdiff(anno.df$name2, expr.df$Gene.Name)
appendgenes.df <- as.data.frame(matrix(data=NA, nrow=length(toaddgenes), 
                                       ncol=ncol(expr.df),
                                       dimnames=list(NULL, colnames(expr.df)))
                                )
appendgenes.df$Gene.Name <- toaddgenes

expr.df <- rbind.data.frame(expr.df, appendgenes.df, stringsAsFactors=F)
rm(appendgenes.df, toaddgenes)

if( any(duplicated(expr.df$Gene.Name)) ){
  stop("Duplicated genes.")
}

# na.rm=T only removes NAs in resulting value column (expression values)
df <- reshape2::melt(data=expr.df, na.rm=F, id.vars=c("Gene.Name"))
rm(expr.df)

df$group <- "ND"
df$group[df$value<0.5] <- "NE" 
df$group[df$value>=0.5 & df$value<=10] <- "LE" 
df$group[df$value>10 & df$value<=1000] <- "ME" 
df$group[df$value>1000] <- "HE" 

if( !identical(is.na(df$value), df$group=="ND") ){
  stop("NA values don't match with ND group")
}

if( any(is.na(df$group)) ){
  stop("Missing group.")
}

if( any(!levels(df$group)%in%names(col.v)) ){
  stop("A level of df$group is not in names of col.v.")
} else {
  df$group <- factor(x=df$group, levels=intersect(names(col.v), unique(df$group)))
}

# Values should not have NAs but if na.rm=F, NAs will be removed and
# there will be a warning.
p <- ggplot(data=df, aes(x=variable)) +
  geom_bar(aes(fill=group), position="fill", na.rm=F) + 
  scale_fill_manual(values=col.v[levels(df$group)]) +
  labs(x="Tissue", y="Fraction", fill="Level",
       title=paste0(out.name, "_fractionrelativetouniquegenesinUCSChg19")) + 
  bgr2 +
  theme(axis.text.x=element_text(size=10, angle=45, hjust=1))
 
ggsave(filename=paste0(out.dir, "/", out.name, "_levelOfExp_barplot.pdf"), 
       width=10, height=10)

# rm(list=ls()); gc()



