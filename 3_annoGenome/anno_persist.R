#Take annotation file and then intersect with the persistent contact matrices

#REMINDERS BEFORE RUNNING
#don't forget to modify the included rows of PERSIST.MX, adjusted for testing only 

#keep column names from annotation file in the output tables
############################################################################################################ 
#Set directories

#functions to outsource
source.dir    <- "/Users/ltamon/ProjectBoard/outsource_R"
#annotation file
annoFile.path <- "/Users/ltamon/Database/ucsc_tables"
#HiC_Human21 persistent contact files
persist.dir   <- "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
#main directory for the task
objective.dir <- "/Users/ltamon/ProjectBoard/HiC_contacts_dim/annoGenome"

#Create directories 

#for outputs, anno~persist files)
#dir.create(paste0(objective.dir, "/out_anno_persist/"))
output.dir <- paste0(objective.dir, "/out_anno_persist")  #check if it exists already

############################################################################################################  
#Packages
##CRAN
#install.packages("foreach") 
library("foreach")
#install.packages("data.table") 
library("data.table") #for fread

##BIOCONDUCTOR
#source("https://bioconductor.org/biocLite.R")
#biocLite("IRanges")
library("IRanges")

############################################################################################################
#Set values

NCPU          <- 1

#length of ij bins
bin.length    <- 40000L

#genome version
genome.ver    <- "hg19"
anno.filename <- "hg19anno" #or "hg38anno"

#2(2MB gap) or "05"(0-5 MB gap), refers to minimum gap accepted to classify a contact, 
#two points should be far enough to filter for contacts within a TAD
gc.v          <- c(2, "05")

chr.v         <- c(1, 2) #c(1:22, "X") 

#set overlap minimum of transcript/annotations with ij bins
#average exon length is 2500 based on this paper Sahakyan and Balasubramanian (2016)
#just to make sure that at least 1 exon length is included in the assigned bin
#https://doi.org/10.1186/s12864-016-2582-9
olap.min      <- 2500L #(bp) or set to default for MAC IRanges findOverlap (1L) 

#column names of contacting bins in persistent contact matrices
i             <- "i"
j             <- "j"

#format of chromosome names in annotation file
#important for subsetting annotation file by chromosome
chr.prefix    <- "chr" #can be Chr# or Chromosome# or Chr X

#column names for txStart and txEnd in annotation file
txStart       <- "txStart"
txEnd         <- "txEnd"

############################################################################################################ 
#Functions

#outsource functions
source(paste0(source.dir, "/UTL_NCPU.R"))
source(paste0(source.dir, "/TrantoRextr/GEN_WhichOverlap.R"))

############################################################################################################ 
#load annotation file 
anno.file <- fread(file=paste0(annoFile.path, "/", anno.filename), header=TRUE, data.table=FALSE)

#check for missing values in annotation file
if( sum(apply(anno.file, MARGIN=c(1,2), function(x) any(is.na(x))))!=0 ){
  stop("Missing values in annotation file.")
}

#lengths of transcripts
txSize <- anno.file[,txEnd]-anno.file[,txEnd]
#just checking if txEnd always greater than txStart as assumed in the annotation file
if( sum(txSize<0) ){
  stop("There are txStart > txEnd.")
} 

foreach(gc=gc.v, .inorder=FALSE) %op% {
  
  #work per chromosome, subset annotation file and load persist matrix per chromosome
  ANNO.PERSIST.MX <- foreach(chr=chr.v, .inorder=FALSE, .combine="rbind") %op% {
    
    #load persist matrix
    load(paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))
    persist.mx   <- PERSIST.MX$hits[1:5000,]  #PERSIST.MX$hits - matrix of ij contacts and tally in each 21 tissue
    persist.ntis <- PERSIST.MX$ntis[1:5000]   #vector of number of tissue the ij contact is present
    #same order as PERSIST.MX$hits
    
    #get the vector of relevant bins for that chromosome, unique(i and j bins), only bins contacting
    #bins sorted to increasing values so the ranges are also increasing
    bins.uniq  <- sort(unique(unlist(persist.mx[, c(i, j)])), decreasing=FALSE)
    #end points of 40-KB bins
    bin.end    <- bins.uniq*bin.length
    #start points of 40-KB bins
    bin.start  <- bin.end-bin.length+1
    
    #subset annotation file per chromosome
    subset.anno <- subset(anno.file, chrom==paste0(chr.prefix, chr))
    
    #assign bins based on txStart-txEnd
    #findoverlaps based on actual ranges of bins
    #Note that WhichOverlap(IRanges) returns indices
    #query   <- subset.anno; annotation (to be classified)
    #subject <- bins.uniq (sorted, increasing); bins (classification)
    
    olap.mx <- WhichOverlap(start.query=subset.anno[,txStart], end.query=subset.anno[,txEnd], space.query=rep("a", nrow(subset.anno)),
                              start.subject=bin.start, end.subject=bin.end, space.subject=rep("a",length(bin.end)),
                              maxgap=-1L, minoverlap=olap.min)

    #annotation assigned to bin
    #note that one annotation can be assigned to multiple bins
    col12 <- cbind(subset.anno[olap.mx[,"query"],], HiC21bin=bins.uniq[olap.mx[,"subject"]])
    
    #select all contact pairs containing the annotation (either in i or j bin)
    #make sure indices are unique 
    col345 <- foreach(itr=as.vector(col12[, "HiC21bin"]), .combine="rbind" ) %op% {
      ind.i        <- which(persist.mx[,i]==itr)
      ind.j        <- which(persist.mx[,j]==itr)
      ind.all      <- unique(c(ind.i, ind.j))  #indices of ij contacts in persist matrix that has the bin of interest
      hits.NA      <- apply(persist.mx[ind.all,c(i, j)], MARGIN=c(1,2), function(x) ifelse(x==itr, x<-NA, x<-x))
      hits.noNA.ij <- na.omit(c(rbind(hits.NA[,i], hits.NA[,j])))
      all.val      <- cbind(HiC21binPartner=hits.noNA.ij, ntis=persist.ntis[ind.all])
      apply(all.val, MARGIN=2, function(x) {paste(x, collapse=",")})
    }
    
    ANNO.PERSIST.MX <- cbind(col12, col345)
    save(ANNO.PERSIST.MX, file=paste0(output.dir, "/", "chr", chr, "_min", gc, "Mb_", genome.ver, "_anno_persist_olapMin", olap.min, ".RData"))
    ANNO.PERSIST.MX 
    
  }
  
  save(ANNO.PERSIST.MX, file=paste0(output.dir, "/", "chrALL_min", gc, "Mb_", genome.ver, "_anno_persist_olapMin", olap.min, ".RData"))
}

#rm(list=ls())
