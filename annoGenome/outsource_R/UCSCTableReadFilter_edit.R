#changed Alex's script so that the table object should be supplied not the path

################################################################################
UCSCTableReadFilter <- function( Table = "lib/ucsc_tables/hg38anno",
                                 Filtering.Scheme = "Human.Nuclear.mRNA (NM)" ){

  print("Reading the specified annotation table...", quote=FALSE)
  #anno <- read.table(RefSeq.Anno.Path, header=TRUE, as.is=TRUE)
  anno <- Table
  
  #--------------------------------------------
  if(Filtering.Scheme == "Human.Nuclear.mRNA"){
    print(paste("Filtering using the '",Filtering.Scheme,"' scheme.",sep=""), quote=F) 
    chrom.name <- as.character(anno[,"chrom"])
    # Selecting only the human nuclear chromosomes:
    anno <- anno[which( chrom.name == "chr1"  |
                        chrom.name == "chr2"  |
                        chrom.name == "chr3"  |
                        chrom.name == "chr4"  |
                        chrom.name == "chr5"  |
                        chrom.name == "chr6"  |
                        chrom.name == "chr7"  |
                        chrom.name == "chr8"  |
                        chrom.name == "chr9"  |
                        chrom.name == "chr10" |
                        chrom.name == "chr11" |
                        chrom.name == "chr12" |
                        chrom.name == "chr13" |
                        chrom.name == "chr14" |
                        chrom.name == "chr15" |
                        chrom.name == "chr16" |
                        chrom.name == "chr17" |
                        chrom.name == "chr18" |
                        chrom.name == "chr19" |
                        chrom.name == "chr20" |
                        chrom.name == "chr21" |
                        chrom.name == "chr22" |
                        chrom.name == "chrX"  |
                        chrom.name == "chrY"    ), ]
    # There are gene names starting from NM_ and NR_ as can be checked by:
    #  unique(sapply(as.character(anno[,"name"]), 
    #         FUN=function(i){unlist(strsplit(i,"_"))[1]},
    #         USE.NAMES=F, simplify=T))
    # For instance, there are 37559 NM_ (mRNA coding) and 9335 NR_ (non-coding)
    # genes in hg19.
    #  length(grep("NR_", anno[,"name"]))
    # 37559 entries are mRNA-coding transcripts in hg19
    anno <- anno[grep("NM_", anno[,"name"]),] 
    
  } else {
    
    Filtering.Scheme == "Human.ncRNA (NR)"
    print(paste("Filtering using the '",Filtering.Scheme,"' scheme.",sep=""), quote=F) 
    chrom.name <- as.character(anno[,"chrom"])
    # Selecting only the human nuclear chromosomes:
    anno <- anno[which( chrom.name == "chr1"  |
                          chrom.name == "chr2"  |
                          chrom.name == "chr3"  |
                          chrom.name == "chr4"  |
                          chrom.name == "chr5"  |
                          chrom.name == "chr6"  |
                          chrom.name == "chr7"  |
                          chrom.name == "chr8"  |
                          chrom.name == "chr9"  |
                          chrom.name == "chr10" |
                          chrom.name == "chr11" |
                          chrom.name == "chr12" |
                          chrom.name == "chr13" |
                          chrom.name == "chr14" |
                          chrom.name == "chr15" |
                          chrom.name == "chr16" |
                          chrom.name == "chr17" |
                          chrom.name == "chr18" |
                          chrom.name == "chr19" |
                          chrom.name == "chr20" |
                          chrom.name == "chr21" |
                          chrom.name == "chr22" |
                          chrom.name == "chrX"  |
                          chrom.name == "chrY"    ), ]
    # There are gene names starting from NM_ and NR_ as can be checked by:
    #  unique(sapply(as.character(anno[,"name"]), 
    #         FUN=function(i){unlist(strsplit(i,"_"))[1]},
    #         USE.NAMES=F, simplify=T))
    # For instance, there are 37559 NM_ (mRNA coding) and 9335 NR_ (non-coding)
    # genes in hg19.
    #  length(grep("NR_", anno[,"name"]))
    # 37559 entries are mRNA-coding transcripts in hg19
    anno <- anno[grep("NR_", anno[,"name"]),] 
    
    #--------------------------------------------
  }
  
  print(paste(length(anno[,1])," records are present.",sep=""), quote=FALSE) 
  return(anno)

}
################################################################################
