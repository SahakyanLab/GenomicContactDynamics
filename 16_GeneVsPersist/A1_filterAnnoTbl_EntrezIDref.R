################################################################################
# Clean the annotation table by:
# a. Get only the relevant chromosomes in UCSC annotation tabe (chr1:22, "X", "Y")
# b. Add unique ID for each row in tables (for future analyses e.g. mapping to 
# persist bins)
# c. Separate tables into reference sequence categories. Accession numbers beginning 
# with NM (coding) and NR (non-coding).
# d. Make HUGOsymbol-EntrezID conversion table for future analyses (e.g. KEGG 
# pathway analysis)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
    os = "Mac"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}

output.dir = annofile.dir
### OTHER SETTINGS #############################################################
annofile.v = "hg19anno" #c("hg19anno", "hg38anno")
chrfilter = FALSE
chr.v = paste("chr", c(1:22, "X", "Y"), sep="")
sepRefSeqCat = FALSE
makeEntrezRefTable = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
source(paste0(lib, "/convertGeneKeyType.R"))
source(paste0(lib, "/RefSeqCat.R")) 
source(paste0(lib, "/UCSCTableReadFilter.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(anno.nme in annofile.v){
  
  if(chrfilter==TRUE){
    
    # Read original annotation table 
    anno.file <- fread(file=paste0(annofile.dir, "/", anno.nme, "_orig_allchr"), 
                       header=TRUE, data.table=FALSE, stringsAsFactors=FALSE) 
    
    # Filter chromosomes
    anno.file <- anno.file[anno.file$chrom%in%chr.v, ]
    
    # Add row number as unique ID for each transcript
    anno.file$uniqueID <- 1:nrow(anno.file)
    write.table(x=anno.file, file=paste0(annofile.dir, "/", anno.nme, "_ALL"),
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    
  }
  
  # Separate into reference sequence categories
  if(sepRefSeqCat==TRUE){
    
    # Read annotation table (chr filtered) 
    anno.file <- fread(file=paste0(annofile.dir, "/", anno.nme, "_ALL"), 
                       header=TRUE, data.table=FALSE, stringsAsFactors=FALSE) 
    
    # Check for RefSeq categories in the annotation file and separate
    RefSeqCat(filePath = paste0(annofile.dir, "/", anno.nme, "_ALL"),
              genomeVersion = gsub(x=anno.nme, pattern="anno", replacement=""),
              accessionCol = "name")
    
    nrow.count <- c()
    for( refseq in c("Human.Nuclear.mRNA", "Human.ncRNA (NR)") ){
      
      ifelse(refseq=="Human.Nuclear.mRNA", lab <- "NM", lab <- "NR")
  
      anno.file.refseq <- UCSCTableReadFilter(Table=anno.file,
                                              Filtering.Scheme=refseq)
      nrow.count[refseq] <- nrow(anno.file.refseq)
      
      write.table(x=anno.file.refseq, file=paste0(annofile.dir, "/", anno.nme, 
                                                  "_", lab),
                  sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    }
    
    if( sum(nrow.count)!=nrow(anno.file) ){
      stop("Reference sequence categories don't add up.")
    }
    
  }
  
  # Make HUGOsymbol-EntrezID reference table
  if(makeEntrezRefTable==TRUE){
   
    # Read annotation table (chr filtered)
    anno.file <- fread(file=paste0(annofile.dir, "/", anno.nme, "_ALL"), 
                       header=TRUE, data.table=FALSE, stringsAsFactors=FALSE) 
    
    genelist <- sort(unique(anno.file$name2)) 
    
    hugo.entrezid <- convertGeneKeyType(genes=genelist, 
                                        origKeyType="SYMBOL",
                                        convertTo="ENTREZID",
                                        org="org.Hs.eg.db", 
                                        convTablePath=NULL,
                                        useKEGGAPI=F,
                                        drop.NA=F)
    
    # Check if HUGO symbols in hugo.entrezid are unique
    if( !identical(genelist, unique(hugo.entrezid$SYMBOL)) ){
      stop("Checkpoint 1.")
    }
    
    write.table(hugo.entrezid, file=paste0(output.dir, "/keyType_conversion/", 
                                           anno.nme, "_SYMBOLtoENTREZID_052020"),
                sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}

# rm(list=ls())

# hg19anno_orig_allchr
# > unique(anno.file$chrom)
#[1] "chr1"                  "chr2"                  "chr3"                  "chr4"                 
#[5] "chr5"                  "chr6"                  "chr7"                  "chr8"                 
#[9] "chr9"                  "chrX"                  "chrY"                  "chr10"                
#[13] "chr11"                 "chr12"                 "chr13"                 "chr14"                
#[17] "chr15"                 "chr16"                 "chr17"                 "chr18"                
#[21] "chr19"                 "chr20"                 "chr21"                 "chr22"                
#[25] "chr6_apd_hap1"         "chr6_cox_hap2"         "chr6_dbb_hap3"         "chr6_mcf_hap5"        
#[29] "chr6_qbl_hap6"         "chr4_ctg9_hap1"        "chr6_mann_hap4"        "chr6_ssto_hap7"       
#[33] "chrUn_gl000211"        "chrUn_gl000212"        "chrUn_gl000213"        "chrUn_gl000215"       
#[37] "chrUn_gl000218"        "chrUn_gl000219"        "chrUn_gl000220"        "chrUn_gl000222"       
#[41] "chrUn_gl000223"        "chrUn_gl000227"        "chrUn_gl000228"        "chr17_ctg5_hap1"      
#[45] "chr1_gl000191_random"  "chr1_gl000192_random"  "chr4_gl000193_random"  "chr4_gl000194_random" 
#[49] "chr7_gl000195_random"  "chr17_gl000205_random" "chr19_gl000209_random"