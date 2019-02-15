################################################################################
#Use GeneVsNtis files to do GO and KEGG enrichment 

################################################################################
lib <- "/Users/ltamon/ProjectBoard/lib"
#lib <- "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
#objective.dir <- "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
#list of genes with HiC strength bins
genelist.dir <- paste0(objective.dir, "/out_GeneVsNtis")
bgrdAll.dir <- paste0(objective.dir, "/out_addEntrezidToAnno")
bgrdContact.dir <- paste0(objective.dir, "/out_backgroundgenes")
#dir.create(paste0(objective.dir, "/C3_funxAnnotation/out_GO"))
output.GO.dir <- paste0(objective.dir, "/B3_funxAnnotation/out_GO")

library(data.table) #for fread
library(foreach)

#package for Installing and Managing Bioconductor Packages
#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("org.Hs.eg.db", version = "3.8")
library(org.Hs.eg.db)

#genome version of annotation file used to derive gene list
genome.ver <- "hg19"

#2(2MB gap) or "05"(0.5 MB gap), refers to minimum gap between contacting bins, 
#two points should be far enough to filter for contacts within a TAD
gc.v <- c("05")

chr.v <- c("ALL") #c(1:22, "X", "ALL") 

#distinguish between protein-coding (NM) and non-coding annotations (NR) or
#NOT (ALL)
data.v <- c("ALL", "NM", "NR")

ntis.v <- c("1", "21")

#Domains/subontologies:
# Biological Process (BP)
# Cellular Component (CC)
# Molecular Function (MF)
ont.v <- c("BP","CC","MF")

#for enrichGO
#organism
OrgDb=org.Hs.eg.db
#type of input genes
keyType="SYMBOL"

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
    foreach (data=c("ALL", "NM", "NR"), .inorder=TRUE) %do% {
      for(background in background.v){
        #set background genes
        if(background=="bgrdAll"){
          uni <- fread(file=paste0(bgrdAll.dir, "/", 
                                   genome.ver, "anno_uniqueHGNC.txt"),
                       data.table=FALSE, header=FALSE)
          uni <- uni[,1]
        }
        if(background=="bgrdContact"){
          uni <- fread(file=paste0(bgrdContact.dir,  "/", genome.ver, "_min", gc, 
                                   "Mb_chr", chr, "_", data, "_background.txt"),
                       data.table=FALSE, header=FALSE)
          uni <- uni[,1]
        }
        if(background=="bgrdPscore1"){
          uni <- fread(file=paste0(genelist.dir,  "/", genome.ver, "_min", gc, 
                                   "Mb_chr", chr, "_", data, "_PScore1", 
                                   "_GeneVsNtis.txt"),
                       data.table=FALSE, header=FALSE)
          uni <- uni[,1]
        }
        #-------------------------------------
        for(i in ntis.v){
          genelist <- fread(file=paste0(genelist.dir, "/", genome.ver, "_min", gc, 
                                        "Mb_chr", chr, "_", data, "_PScore", i,  
                                        "_GeneVsNtis.txt"),
                            data.table=FALSE, header=FALSE)
          genelist <- as.character(genelist[,1])
          #--------------------
          #GO enrichment analysis   
          prefix <- paste0(output.GO.dir, "/GO_", genome.ver, "_min", gc, 
                           "Mb_chr", chr, "_", data, 
                           "_PScore", i, "_", background)
          
          if(length(ont.v)==1){
            suffix <- paste0("_", ont.v)
          } else {
            suffix <- paste0("_", paste(ont.v, collapse="_"))
          }
          pdf(height=8, width=10, file=paste0(prefix, suffix, ".pdf"))
          for(ont in ont.v){
            ego <- enrichGO(gene=genelist,
                            OrgDb=OrgDb,
                            keyType=keyType,
                            ont=ont,
                            pvalueCutoff=0.05,
                            pAdjustMethod="BH",
                            universe=uni)
            write.table(as.data.frame(ego), file=paste0(prefix, "_", ont, ".txt"),
                        col.names=TRUE, row.names=FALSE, quote=FALSE)
            print(dotplot(ego, x="count", showCategory=30))
          }  #ont for loop
          dev.off() 
        }  #ntis for loop
      } #background for loop
    } #data foreach end loop
  }
}

#rm(list=ls())
 





