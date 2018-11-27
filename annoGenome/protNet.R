#Take annotations(HUGOnames)~ntis(ntis for number of tissues, referring to
#presence of persistent contacts in 21 tissues) gene ontology and pathway analysis

############################################################################################################ 
#Set directories

#main directory for the objective or project
objective.dir <- "/Users/ltamon/ProjectBoard/HiC_contacts_dim/annoGenome"
#directory of HUGOnames~ntis data
geneList.dir  <- paste0(objective.dir, "/out_HUGOname_ntis") 

############################################################################################################ 
#Packages
##CRAN 

#package for Installing and Managing Bioconductor Packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("STRINGdb", version = "3.8")
library(STRINGdb)

############################################################################################################ 
#Set values
#genome version of annotation file used to derive gene list
genome.ver <- "hg19"

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gc.v       <- c(2, "05")

#name of column containing the annotations, HUGO names in this case (refer to gene list file)
GENE       <- paste0(genome.ver,"HUGOname")

#name of column containing number of tissues of persistent contacts
#for each annotation (refer to gene list file)
NTIS       <- "ntis"

############################################################################################################
#Functions

#Protein network
#using STRINGdb (database of known and predicted protein-protein interac- tions)
#The interactions include direct (physical) and indirect (functional) associations

get.protnet <- function(genes = genes,
                        plotname = plotname){
  
  df <- data.frame(genes=genes)
  #map the genes to the STRING database identifiers
  #$map(my_data_frame, my_data_frame_id_col_names, takeFirst = , removeUnmappedRows = , quiet = )
  #quiet=FALSE, issue warning about unmapped values
  str.mapped <- string_db$map( df, "genes", removeUnmappedRows=TRUE, quiet=FALSE )
  pdf(height=8, width=8, file=plotname)
  #extract genes mapped (specify how many) and produce an image of the STRING network for those
  string_db$plot_network(
    str.mapped$STRING_id[1:min(400, length(str.mapped$STRING_id))]
  )
  dev.off()
}

#-----------------------------------------------
#for testing
gc = 2
i=1

#load geneList differentiating between NM and NR
#loaded object is HUGO.NTIS.NMNR
load( paste0(geneList.dir, "/chrALL_min", gc, "Mb_", genome.ver, "_NMNRHUGO_ntis.RData") )
HUGO.NTIS.NM <- subset(HUGO.NTIS.NMNR, HUGO.NTIS.NMNR[,"Grp"]=="NM")
ntis.uniq    <- as.character(unique(HUGO.NTIS.NM[,NTIS]))

string_db    <- STRINGdb$new(version="10", species=9606,
                             score_threshold=400, input_directory="")
#score_threshold
#a threshold for the combined scores of the interactions, such that 
#any interaction below that threshold is not loaded in the object (by default 
#the score threshold is set to 400!)

#subset per ntis
for(i in ntis.uniq){
  sbset <- subset(HUGO.NTIS.NM, HUGO.NTIS.NM[,NTIS]==i)
  get.protnet(genes=as.character(sbset[,GENE]), 
              plotname=paste0(geneList.dir, "/min", gc, "Mb_", genome.ver, "_", i, "tis_protNet.pdf"))
}
#rm(list=ls())




