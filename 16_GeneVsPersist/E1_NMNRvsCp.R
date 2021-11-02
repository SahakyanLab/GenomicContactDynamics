################################################################################
# Fraction of NM (coding genes) and NR (non-coding genes) vs. Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/16_GeneVsPersist"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
annofilePath = paste0(annofile.dir, "/hg19anno_ALL")
# Use longest transcript to remove bias from uneven count of transcripts per gene
#annofilePath = paste0(annofile.dir, "/hg19annoLTr_ALL")
# Unique HUGO genes overlapping with Cp regions (any type of overlap)
CpGenesPath = paste0(wk.dir, "/out_anno_union/min2Mb_ALL_name2")
out.dir = paste0(wk.dir, "/out_NMNRvsCp")
### OTHER SETTINGS #############################################################
# Annotation file prefix
anno.nme = "hg19anno" #"hg19annoLTr"
gcb="min2Mb"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
source(paste0(paste0(lib, "/GG_bgr.R")))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
anno.df <- read.delim(file=annofilePath, header=TRUE, stringsAsFactors=FALSE)
# Non coding genes (NR_)
NRgenes <- anno.df[anno.df$cdsStart==anno.df$cdsEnd, "name2"]

# Retrieve Cp genes
temp <- readLines(con=CpGenesPath)
CpGenesInd <- grep(x=temp, pattern=">all_genes_cp_", fixed=TRUE)
CpGenes <- temp[CpGenesInd + 1L]
names(CpGenes) <- gsub(x=temp[CpGenesInd], pattern=">all_genes_cp_|_end", replacement="")
rm(temp, CpGenesInd); gc()
CpGenes.len <- length(CpGenes)
# List of genes overlapping with all long-range contact regions and regions
# of Cp 1 to 21. 
CpGenes <- strsplit(x=CpGenes, split=";", fixed=TRUE)

NMNR.DF <- sapply(X=names(CpGenes), simplify=FALSE, FUN=function(cp){
  if( any(duplicated(CpGenes[[cp]])) ){ stop(paste0("Duplicated genes in ", cp)) }
  n0 <- length(CpGenes[[cp]])
  frNR <- round(sum(CpGenes[[cp]]%in%NRgenes)/n0, digits=3)
  frNM <- round(sum(!CpGenes[[cp]]%in%NRgenes)/n0, digits=3)
  cp <- as.character(cp)
  if( (frNR+frNM)==1 ){
    print(paste0(cp, " done!"), quote=FALSE)
    return( cbind.data.frame(Cp=c(cp,cp), type=c("NM", "NR"), fr=c(frNM, frNR)) )
  } else {
    stop("Sum is not 1.")
  }
})
NMNR.DF <- do.call("rbind.data.frame", NMNR.DF)

# Convert Cp column to factor, arranging the levels for the plot
NMNR.DF$Cp <- as.character(NMNR.DF$Cp)
NMNR.DF$Cp <- factor(x=NMNR.DF$Cp, levels=unique(NMNR.DF$Cp))

ggplot(data=NMNR.DF, aes(x=Cp, y=fr, fill=type)) +
  geom_col() +
  labs(title=paste0(gcb, "_", anno.nme), x=expression(bold("c"["p"])),
       y="Fraction", fill=NULL) + 
  bgr2 +
  theme(axis.text.x = element_text(size=10, angle=360, colour="black"))
ggsave(filename=paste0(out.dir, "/", gcb, "_", anno.nme, "_NMNRvsCp.pdf"),
       width=10, height=10)

# rm(list=ls()); gc()

