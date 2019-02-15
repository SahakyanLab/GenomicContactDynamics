#Use HGNCsymbVsPersist files to do GO and KEGG enrichment 

################################################################################
#Set directories

#functions to source
lib <- "/Users/ltamon/ProjectBoard/lib"
#lib <- "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
#objective.dir <- "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
#list of genes with HiC strength bins
genelist.dir <- paste0(objective.dir, "/out_GeneVsNtis")
bgrdAll.dir <- paste0(objective.dir, "/out_addEntrezidToAnno")
bgrdContact.dir <- paste0(objective.dir, "/out_backgroundgenes")
#dir.create(paste0(objective.dir, "/C3_funxAnnotation/out_KEGG"))
output.KEGG.dir <- paste0(objective.dir, "/B3_funxAnnotation/out_KEGG")

################################################################################
#Packages
##CRAN 
#install.packages("data.table") 
library(data.table) #for fread
library(foreach)

#package for Installing and Managing Bioconductor Packages
#install.packages("BiocManager")

#BiocManager::install("clusterProfiler", version="3.10.1")
library(clusterProfiler)
#BiocManager::install("KEGG.db")

################################################################################
#Set values

#genome version of annotation file used to derive gene list
genome.ver <- "hg19"

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gc.v <- c(2, "05")

chr.v <- c("ALL") #c(1:22, "X", "ALL") 

#distinguish between protein-coding (NM) and non-coding annotations (NR) or
#NOT (ALL)
data.v <- c("ALL", "NM", "NR")

ntis.v <- 1:21

#name of column containing the annotations, 
#HGNC symbols in this case (refer to gene list file)
gene <- "Gene"

#for enrichKEGG
#organism
org="hsa"
#type of input genes
keyType="ncbi-geneid"
#library(KEGG.db)

#set background for analysis
background.v <- c("bgrdAll", "bgrdContact", "bgrdPscore1")
#allgenes <- take all annotations(HGNC symbols) from annotation file
#bgrdContact <- all annotations associated with contacts, distinguishing
#between data types ("ALL", "NM", "NR")
#bgrdPscore1 <- Persistence score 1 genes

#gc=2
#chr="ALL"
#data="ALL"
#i="21"
#background <- "bgrdContact"

################################################################################
################################################################################
foreach(gc=gc.v, .inorder=TRUE) %do% {
  foreach(chr=chr.v, .inorder=TRUE) %do% {
    foreach(data=c("ALL", "NM", "NR"), .inorder=TRUE) %do% {
      
      for(background in background.v){
        #set background genes
        
        if(background=="bgrdAll"){
          uni <- fread(file=paste0(bgrdAll.dir, "/", 
                                   genome.ver, "anno_uniqueHGNC_entrezId.txt"),
                       data.table=FALSE, header=TRUE)
          uni <- uni[,"entrezgene"]
          uni <- na.omit(uni)
        }
        
        if(background=="bgrdContact"){
          uni <- fread(file=paste0(bgrdContact.dir,  "/entrezId_", genome.ver, "_min",
                                   gc, "Mb_chr", chr, "_", data, "_background.txt"),
                       data.table=FALSE, header=TRUE)
          uni <- uni[,"entrezgene"]
          uni <- na.omit(uni)
        }
        
        if(background=="bgrdPscore1"){
          uni <- fread(file=paste0(bgrdContact.dir,  "/entrezId_", genome.ver, "_min",
                                   gc, "Mb_chr", chr, "_", data, "_background.txt"),
                       data.table=FALSE, header=TRUE)
          uni <- uni[,"entrezgene"]
          uni <- na.omit(uni)
        }
        
        #-------------------------------------
        for(i in ntis.v){
          
          genelist <- fread(file=paste0(genelist.dir, "/entrezId_", genome.ver, "_min", 
                                        gc, "Mb_chr", chr, "_", data, "_PScore", i,  
                                        "_GeneVsNtis.txt"),
                            data.table=FALSE, header=TRUE)
          genelist <- na.omit( as.character(genelist[,"entrezgene"]) )
          
         
          #--------------------
          #KEGG enrichment analysis   
          
          ekegg <- enrichKEGG(gene=genelist,
                              organism=org,
                              #one of "kegg", 'ncbi-geneid', 'ncib-proteinid' 
                              #and 'uniprot'
                              keyType=keyType, #NCBI Entrez ID
                              pAdjustMethod = "BH",
                              universe=uni,
                              #defaults
                              pvalueCutoff=0.05, minGSSize=10, 
                              maxGSSize=500, qvalueCutoff=0.2,
                              #FALSE: download the latest KEGG data 
                              #for enrichment analysis
                              #TRUE: use KEGG.db
                              use_internal_data=FALSE         
          )
          
          write.table(as.data.frame(ekegg), 
                      file=paste0(output.KEGG.dir, "/KEGG_", genome.ver, "_min", 
                                  gc, "Mb_chr", chr, "_", data, "_PScore",
                                  i, "_", background, ".txt"),
                      col.names=TRUE, row.names=FALSE, quote=FALSE)
          
          #viewKEGG(ekegg, pathwayID, foldChange, color.low = "green",
          #         color.high = "red", kegg.native = TRUE,
          #         out.suffix = ".pdf")
          
        } #ntis for loop
      } #background for loop
    }
  }
}

#rm(list=ls())
 
 






