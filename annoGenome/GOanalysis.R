#Take annotations(HUGOnames)~ntis(ntis for number of tissues, referring to
#presence of persistent contacts in 21 tissues) gene ontology and pathway analysis

############################################################################################################ 
#Set directories

#main directory for the task
objective.dir <- "/Users/ltamon/ProjectBoard/HiC_contacts_dim/annoGenome"
#annotation file
annoFile.path <- "/Users/ltamon/Database/ucsc_tables"
#HUGOnames~ntis data
geneList.dir  <- paste0(objective.dir, "/out_HUGOname_ntis") 

############################################################################################################ 
#Packages
##CRAN 

#install.packages("data.table") 
library("data.table") #for fread

#package for Installing and Managing Bioconductor Packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("AnnotationDbi", version = "3.8")
library(AnnotationDbi)

#BiocManager::install("org.Hs.eg.db", version = "3.8")
library(org.Hs.eg.db)
#citation("org.Hs.eg.db")
# Marc Carlson (2018). org.Hs.eg.db: Genome wide annotation for
# Human. R package version 3.7.0.

#BiocManager::install("clusterProfiler", version = "3.8")
library(clusterProfiler) #for enrichGO function
#clusterProfiler v3.8.1  For help: https://guangchuangyu.github.io/software/clusterProfiler
#citation
# Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. 
# clusterProfiler: an R package for comparing biological themes among gene clusters. 
# OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.

############################################################################################################ 
#Set values
#genome version of annotation file used to derive gene list
genome.ver    <- "hg19"
anno.filename <- "hg19anno" #or "hg38anno"

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gc.v    <- c(2, "05")

#name of column containing the annotations, HUGO names in this case (refer to gene list file)
GENE    <- paste0(genome.ver,"HUGOname")

#name of column containing number of tissues of persistent contacts
#for each annotation (refer to gene list file)
NTIS    <- "ntis"

#column name of HUGO names in annotation file (for specifiying background for GO analysis)
name2   <- "name2"

#organism, specify OrgDb object
OrgDb   <- org.Hs.eg.db

#keytype of input gene and make sure it's provided by OrgDb
#keytypes(org.Hs.eg.db) #function from AnnotationDbi package 
keyType <- "SYMBOL" #(e.g. "ENTREZID", "ENSEMBL")

#set background for GO analysis
#take only annotations present in at least 1 of the tissues (ntis>=1)
ANNOpresent <- FALSE 
#if  FALSE, take all annotations(HUGO symbols) from annotation file

############################################################################################################
#Functions

#GO Analysis
#enrichGO uses VECTOR of genes

enGO <- function(genes = geneList,
                 filename = filename,
                 bgr = uni){

    # Biological Process (BP)
    # Cellular Component (CC)
    # Molecular Function (MF)

    pdf(height=8, width=10, file=filename)
    for(ont in c("BP","CC","MF")){
        ego <- enrichGO(gene=genes,
          OrgDb=OrgDb,
          keyType=keyType,
          ont=ont,
          pvalueCutoff=0.05,
          pAdjustMethod="BH",
          universe=bgr)
          print(dotplot(ego, x="count", showCategory=30))
    }
    dev.off()
}

#-----------------------------------------------
gc=2
i=1

#load gene list (HUGO.NTIS)
load( paste0(geneList.dir, "/chrALL_min", gc, "Mb_", genome.ver, "_ALLHUGO_ntis.RData") )
ntis.uniq <- as.character(unique(HUGO.NTIS[,NTIS])) 

if (ANNOpresent==TRUE){
  uni <- as.character(unique(HUGO.NTIS[,GENE]))
} else {
  #load annotation file to specify background
  anno.file <- fread(file=paste0(annoFile.path, "/", anno.filename), header=TRUE, data.table=FALSE)
  uni       <- as.character(unique(anno.file[,name2]))
}

#subset per ntis
for(i in ntis.uniq) {
  sbset <- subset(HUGO.NTIS, HUGO.NTIS[,NTIS]==i)
  enGO(genes=as.character(sbset[,GENE]), bgr=uni,
       filename=paste0(geneList.dir, "/min", gc, "Mb_", genome.ver, "_", i, "tis_enGO_BP_CC_MF.pdf"))
}

#rm(list=ls())




