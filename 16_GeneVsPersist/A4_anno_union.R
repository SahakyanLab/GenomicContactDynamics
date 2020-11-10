################################################################################
# Union of genes from all chr
# Note: Header with no genes still have an additional line containing "" 
# not NULL or NA, length==1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
    annofile.dir = "/t1-data/user/ltamon/Database/ucsc_tables/hsa_geneAnno"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(wk.dir, "/out_anno_perChr")
out.dir = paste0(wk.dir, "/out_anno_union")
# Converter of HUGO gene symbols to ncbi-geneid for KEGG
hugoEntrezPath = paste0(annofile.dir, "/keyType_conversion/hg19anno_SYMBOLtoENTREZID_052020")
### OTHER SETTINGS #############################################################
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
anno.nme = "hg19anno"
suffix = "name2" # "name2" | "uniqueID"
gcb = "min05Mb" # "min2Mb" | "min05Mb"
refseq.v = c("ALL", "NM", "NR")
chr.v = paste("chr", c(1:22, "X"), sep="")
nCPU = 2L
### FUNCTION ###################################################################
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools) # isplitVector
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, " ", suffix), quote=FALSE)

chr.v.len <- length(chr.v)

if(suffix=="name2"){
  # HUGOsymbol-entrezID conversion table
  tbl <- fread(file=hugoEntrezPath,
               data.table=FALSE, header=TRUE, stringsAsFactors=FALSE)
}

for(refseq in refseq.v){
  
  if(suffix=="name2"){
    # All unique genes from annotation table (chr filtered)
    # To be added later to the final text file
    allgenes <- unique(fread(file=paste0(annofile.dir, "/", anno.nme, "_", refseq), 
                             header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)[["name2"]])
  }
  
  compareToGenes <- list()
  even.lnes <- NULL
  even.lnes.len <- NULL
  
  for(c in 1:chr.v.len){
    chr <- chr.v[c]
    fle <- paste0(data.dir, "/", chr, "_", gcb, "_", refseq, "_", suffix)
    lneslist <- readLines(con=fle) 
    # Change blank lines to "none"
    lneslist[lneslist==""] <- "none"
    lneslist <- strsplit(x=lneslist, split=";")
    
    if( c!=1 ){
      
      #### FOREACH EXECUTION #########
      out.lst <- foreach(itr=isplitVector(x=1:even.lnes.len, chunks=nCPU),
                         .inorder=TRUE, .combine="c",
                         .export=c("even.lnes", "compareToGenes", "lneslist"), 
                         .noexport=ls()[!ls()%in%c("even.lnes", "compareToGenes", "lneslist")]
      ) %op% {
        
        chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
          # Union (should remove duplicates)
           genes.uniq <- union(compareToGenes[[ even.lnes[i] ]], 
                               lneslist[[ even.lnes[i] ]])
          # Unique
          #genes.uniq <- unique( c(compareToGenes[[ even.lnes[i] ]],
          #                        lneslist[[ even.lnes[i] ]]) )
        })
        return(chunk)
      } # foreach loop end
      
      compareToGenes[even.lnes] <- out.lst
        
      ### END OF FOREACH EXECUTION ###
      
    } else {
      # Only lines with HUGO/uniqueID
      lneslist.len <- length(lneslist)
      even.lnes <- seq(from=2, to=lneslist.len, by=2)
      even.lnes.len <- length(even.lnes)
      compareToGenes <- lneslist
    }
    
    rm("lneslist"); gc()
    
    print(chr)
    
  } # chr end for loop
  
  compareToGenes <- lapply(X=compareToGenes, FUN=function(x){
    # Collapse genes to string
    # Remove "none" in gene strings but keep blank lines (no genes at all)
    # as "none"
    gsub(x=paste(x, collapse=";"), pattern="none;|;none", replacement="")
  })
  # Transform to vector 
  compareToGenes <- unlist(compareToGenes, use.names=FALSE)
  
  if(suffix=="name2"){
    # All genes that have non-0 counts in at least one of 21 Hi-C datasets 
    CpSbset <- compareToGenes[grep(x=compareToGenes, pattern="all_genes_cp_", fixed=TRUE)+1L]
    #CpSbset <- compareToGenes[seq(from=2, to=42, by=2)]
    # Remove "none" lines (no genes at all)
    CpSbset <- CpSbset[CpSbset!="none"]
    allgenesCp <- unique( 
      unlist( strsplit(x=CpSbset, split=";"), use.names=FALSE )
    )
    
    # Append 
    compareToGenes <- append( x=compareToGenes, 
                              values=c(">all_genes_end", 
                                       paste(allgenes, collapse=";"), 
                                       ">all_genes_cp_HiC_all_end", 
                                       paste(allgenesCp, collapse=";")),
                              after=0 )
    
    rm("CpSbset", "allgenesCp")
  }
  
  write(x=compareToGenes, file=paste0(out.dir, "/", gcb, "_", refseq, "_", suffix),
        append=TRUE)

  if(suffix=="name2"){
    
    # EntrezID version of txt files
    # Convert HUGO symbols to entrezID
    
    # Indices of gene strings
    ind <- grep(pattern=">all_genes", x=compareToGenes, fixed=TRUE) + 1L
    compareTolist <- strsplit(x=compareToGenes[ind], split=";")
    lst <- lapply(X=compareTolist, FUN=function(x){
      x <- unique(tbl[tbl$SYMBOL%in%x, "ENTREZID"])
      paste(x[!is.na(x)], collapse=";")
    })
    compareToGenes[ind] <- unlist(x=lst, use.names=FALSE)
    
    rm("ind", "compareTolist", "lst"); gc()
    
    write(x=compareToGenes, file=paste0(out.dir, "/", gcb, "_", refseq, "_entrezID"),
          append=TRUE)
  }
  
  rm("compareToGenes"); gc()

} # refseq.v for loop end

# rm(list=ls())




