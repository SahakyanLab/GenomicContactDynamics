################################################################################
# Do some checks of the final mutation datasets. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
chrLenFile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
src.dir = paste0(wk.dir, "/out_filter")
out.dir = paste0(wk.dir, "/out_check")
### OTHER SETTINGS #############################################################
src.id = "hg38ToHg19" # "hg38ToHg19" | "Hg19"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
source(paste0(wk.dir, "/lib/exploreData.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(file=paste0(src.dir, "/CosmicNCV_", src.id, "_final.RData"))

x <- list()

# Check if mutations are within GRCh37/hg19 chromosome length
if( any(c(ncv.df$start, ncv.df$end)<1) ){
  stop("Zero/Negative start/end coordinate.")
}
chrLen.df <- read.delim(file=chrLenFile, header=TRUE, stringsAsFactors=FALSE)
chr.v <- na.omit(unique(ncv.df$chr))
outChr.TF <- unlist(sapply(X=chr.v, simplify=TRUE, FUN=function(chr){
  chr.len <- chrLen.df[chrLen.df$chromosome==chr,"length.bp"]
  return( max(ncv.df[ncv.df$chr==chr,"end"], na.rm=TRUE)>chr.len )
}))
if( any(unname(outChr.TF)) ){
  stop("Mutations outside hg19 chromosome length.")
}
rm(chrLen.df, outChr.TF, chr.v)

  # Sample ID and name
ulen1 <- length(unique(ncv.df$ID_SAMPLE))
ulen2 <- length(unique(ncv.df$`Sample name`))
ulen3 <- length(unique(paste0(ncv.df$`Sample name`, ncv.df$ID_SAMPLE)))
if( length(unique(c(ulen1, ulen2, ulen3)))!=1 ){
  stop("Problem with samples.")
}

# SBS mutation type
if( any(!ncv.df$MUT%in%c("T>G","T>C","T>A","C>A","C>G","C>T")) ){
  stop("Problematic WT_SEQ/MUT_SEQ.")
}

# GENOMIC_MUTATION_ID
if( any(duplicated(na.omit(ncv.df[,c("Sample name", "ID_SAMPLE", "GENOMIC_MUTATION_ID")]))) ){
  x$dupGMIPerSample <- TRUE
} else {
  x$dupGMIPerSample <- FALSE
}

# No sample has >1 mutations in a site
if( any(duplicated(na.omit(ncv.df[,c("Sample name", "ID_SAMPLE", "genome position")]))) ){
  x$multMutPerPosPerSample <- TRUE
} else {
  x$multMutPerPosPerSample <- FALSE
}

# LEGACY_MUTATION_ID
if( any(duplicated(na.omit(ncv.df[,c("Sample name", "ID_SAMPLE", "LEGACY_MUTATION_ID")]))) ){
  x$dupLMIPerSample <- TRUE
} else {
  x$dupLMIPerSample <- FALSE
}

# Match of added chr, start and end coordinates to genome position
position <- strsplit(x=ncv.df$`genome position`, split=":|-")

for( y in c("chr", "start", "end") ){
  
  v <- unlist(lapply(X=position, FUN=function(x)x[ c(chr=1, start=2, end=3)[y] ] ))
  
  if(y=="chr"){
    v <- paste0("chr", v)
    v[v=="chr23"] <- "chrX"
    v[v=="chr24"] <- "chrY"
    v[v=="chr25"] <- "chrMT"
  }
  
  x[[y]] <- c(same=NA, diff=NA)
  x[[y]]["same"] <- sum(v==ncv.df[[y]])
  x[[y]]["diff"] <- sum(v!=ncv.df[[y]])
  
}

# Only coordinate that don't match in hg19 is because the start
# and end were written in exponential. See below.

#-------------------
# Change genome position 
ncv.df$`genome position` <- paste0(ncv.df$chr, ":", ncv.df$start, "-", ncv.df$end)
ncv.df$`genome position` <- gsub(x=ncv.df$`genome position`, pattern="chr", 
                                 replacement="", fixed=TRUE)
# Be aware of mutations with missing coordinates after liftover
ncv.df[is.na(ncv.df$start) & is.na(ncv.df$end),] <- NA
#ncv.df[1:2,] <- NA # REMOVE

x <- c(x, exploreData(df=ncv.df))

x$MUT <- table(ncv.df$MUT, useNA="always")/x$totmut
x$start <- x$start/x$totmut
x$end <- x$end/x$totmut
x$Coding <- table(ncv.df$Coding, useNA="always")/x$totmut

save(x, file=paste0(out.dir, "/CosmicNCV_", src.id, "_check.RData"))

# rm(list=ls()); gc()

#ncv.df[v!=ncv.df[[y]],]
#Sample name ID_SAMPLE ID_tumour                       Primary site
#1765212 006-0018-01TD   1921501   1808819 haematopoietic_and_lymphoid_tissue
#Site subtype 1 Site subtype 2 Site subtype 3 Primary histology
#1765212             NS             NS             NS lymphoid_neoplasm
#Histology subtype 1
#1765212 chronic_lymphocytic_leukaemia-small_lymphocytic_lymphoma
#Histology subtype 2 Histology subtype 3 GENOMIC_MUTATION_ID
#1765212                  NS                  NS        COSV65186327
#LEGACY_MUTATION_ID zygosity GRCh       genome position
#1765212       COSN18983707  Unknown   37 3:155000000-155000000
#Mutation somatic status WT_SEQ MUT_SEQ SNP
#1765212 Confirmed somatic variant      G       T   n
#FATHMM_MKL_NON_CODING_SCORE FATHMM_MKL_NON_CODING_GROUPS
#1765212                     0.15257                           NA
#FATHMM_MKL_CODING_SCORE FATHMM_MKL_CODING_GROUPS Whole_Genome_Reseq
#1765212                 0.01533                       NA                  y
#Whole_Exome ID_STUDY PUBMED_PMID            HGVSG      end    start
#1765212           n      340    21642962 3:g.155000000G>T 1.55e+08 1.55e+08
#chr MUT
#1765212 chr3 C>A
