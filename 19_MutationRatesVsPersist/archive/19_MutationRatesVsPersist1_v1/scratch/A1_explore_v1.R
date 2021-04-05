################################################################################
# Explore COSMIC non-coding mutation data.
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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
genome = 38
CosmicNCVPath = paste0(data.dir, "/cosmic/GRCh", genome, "/CosmicNCV.GRCh", 
                       genome, ".tsv")
out.dir = paste0(wk.dir, "/out_explore")
### OTHER SETTINGS #############################################################
sample.len = 10000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
x <- list()

ncv.df <- fread(file=CosmicNCVPath, header=TRUE, data.table=FALSE) 

# Total mutations
x$total <- length(ncv.df[,1])

x$columns <- colnames(ncv.df)

# Check for NAs 
tmp <- list()
for( i in 1:ncol(ncv.df) ){
  colnme <- colnames(ncv.df)[i]
  if( any(is.na(ncv.df[[colnme]])) ){
    tmp[[colnme]] <- colnme
  }
  rm(colnme)
}
x$withNAs <- names(tmp); rm(tmp)

# Primary site
x$`Primary site` <- table(ncv.df$`Primary site`)

# Primary histology
x$`Primary histology` <- table(ncv.df$`Primary histology`)

# Zygosity
x$zygosity <- table(ncv.df$zygosity)

# Genome reference
x$GRCh <- table(ncv.df$GRCh)

# Genomic position
tmp <- strsplit(x=ncv.df$`genome position`, split=":|-")
tmp1 <- unlist(lapply( X=1:length(tmp), FUN=function(i){
  x <- as.numeric(tmp[[i]])
  (x[3]-x[2]+1L)==nchar(ncv.df$WT_SEQ[i])
}))
x$onebased <- table(tmp1)
rm(tmp1)

# Chromosome
x$chr <- paste0("chr", 
                unique( unlist(lapply(X=tmp, FUN=function(x) x[1])) )
                )
rm(tmp)

# Mutation somatic status
x$`Mutation somatic status` <- table(ncv.df$`Mutation somatic status`)

# WT SEQ --> MUT SEQ
x$wtmutdiff <- sum(ncv.df$WT_SEQ!=ncv.df$MUT_SEQ)
x$wtmutdifflen <- sum( nchar(ncv.df$WT_SEQ)!=nchar(ncv.df$MUT_SEQ) )
x$WT_SEQ <- table(ncv.df$WT_SEQ)
x$MUT_SEQ <- table(ncv.df$MUT_SEQ)

# SNP
x$SNP <- table(ncv.df$SNP)

# FATHMM_MKL_NON_CODING_SCORE & FATHMM_MKL_CODING_SCORE
for( y in c("FATHMM_MKL_NON_CODING_SCORE", "FATHMM_MKL_CODING_SCORE") ){
  
  tmp <- ncv.df[[y]]
  tmp <- tmp[!is.na(tmp)]
  x[[y]] <- c(Missing=sum(is.na(ncv.df[[y]])),
              Neutral=sum(tmp<=0.5),
              Pathogenic=sum(tmp>=0.7), 
              Mid=sum(tmp>0.5 & tmp<0.7))
                                     
  if(sum(x[[y]])!=x$total){
    stop("Checkpoint 1.")
  }
  rm(tmp)
  
}

# FATHMM_MKL_NON_CODING_GROUPS & FATHMM_MKL_CODING_GROUPS
for( y in c("FATHMM_MKL_NON_CODING_GROUPS", "FATHMM_MKL_CODING_GROUPS") ){
 
  x[[y]] <- c( Missing=sum(is.na(ncv.df[[y]])), table(ncv.df[[y]]) )
  
  if(sum(x[[y]])!=x$total){
    stop("Checkpoint 2.")
  }
   
}

# Whole Genome Reseq
x$Whole_Genome_Reseq <- table(ncv.df$Whole_Genome_Reseq)

# Whole_Exome
x$Whole_Exome <- table(ncv.df$Whole_Exome)

#-------------------Remove samples with multiple mutations at the same site
# Custom id for a sample; Need alternative to ID_SAMPLE because there can be 
# multiple ids, if the same sample has been entered into the database multiple 
# times from different papers.

# Also, there can be multiple sample ids linked to the same sample name if it 
# is unclear whether the same sample name between different publications is 
# indeed the same sample.

# TO DO: Remove ambiguity. Take only unique pairs of ID_SAMPLE-Sample name

custom.id <- paste(ncv.df$`Primary site`, ncv.df$`Site subtype 1`, 
                   ncv.df$`Site subtype 2`, ncv.df$`Site subtype 3`,
                   ncv.df$`Primary histology`, ncv.df$`Histology subtype 1`,
                   ncv.df$`Histology subtype 2`, ncv.df$`Histology subtype 3`)

# Trim down dataset
#ncv.df <- ncv.df[, c("Sample name", "ID_SAMPLE", "GENOMIC_MUTATION_ID", 
#                     "LEGACY_MUTATION_ID", "genome position")]

x$multMutInPos <- c(pos_ID_SAMPLE=NA, pos_Sample_name=NA, pos_custom.id=NA, 
                    anyOfSampleIDNAME=NA, anyOf3=NA)

tmp <- paste0(custom.id, ncv.df$`genome position`)
TF1 <- duplicated(tmp) | duplicated(tmp, fromLast=TRUE)
x$multMutInPos["pos_custom.id"] <- sum(TF1)
rm(tmp, custom.id); gc()

tmp <- paste0(ncv.df$ID_SAMPLE, ncv.df$`genome position`)
TF2 <- duplicated(tmp) | duplicated(tmp, fromLast=TRUE)
x$multMutInPos["pos_ID_SAMPLE"] <- sum(TF2)
rm(tmp)

tmp <- paste0(ncv.df$`Sample name`, ncv.df$`genome position`)
TF3 <- duplicated(tmp) | duplicated(tmp, fromLast=TRUE)
x$multMutInPos["pos_Sample_name"] <- sum(TF3)
rm(tmp); gc()

SampDrop.TF <- TF1 | TF2 | TF3

x$multMutInPos["anyOf3"] <- sum(SampDrop.TF)
x$multMutInPos["anyOfSampleIDNAME"] <- sum(TF2 | TF3)
rm(TF1, TF2, TF3); gc()

# COSMIC Identifiers, GENOMIC_MUTATION_ID and LEGACY_MUTATION_ID
# https://cosmic-blog.sanger.ac.uk/searchable-cosmic-identifiers/#:~:text=The%20genomic%20mutation%20identifier%20(COSV,assemblies%20(GRCh37%20and%20GRCh38).

# GENOMIC_MUTATION_ID 
# The genomic mutation identifier (COSV) indicates the definitive position of 
# the variant on the genome. This identifier is trackable and stable between 
# different versions of the release and remains the same between different 
# assemblies (GRCh37 and GRCh38). The genomic mutation identifier is now the 
# preferred way to identify mutations. Note that mutations with no known genomic
# coordinates will not have a value for this identifier, and may be identified 
# using the legacy mutation identifier. The genomic mutation identifier is 
# referred to as GENOMIC_MUTATION_ID in the download files.

# Although it dictates the definitive position of the variant in the genome (independent of
# sample), the GENOMIC_MUTATION_ID for a mutation can have different coordinates depending
# on the transcript. I just don't know whether this is the case for the dataset that we
# have cause I assume that the coordinates used are based on the reference. 

# TO DO: After removing mutations occuring in same site in same sample (based on ID_SAMPLE
# and Sample name), remove mutations with the same GENOMIC_MUTATION_ID per sample.

# LEGACY_MUTATION_ID

# The existing COSM and COSN mutation identifiers are now referred to as legacy mutation 
# identifiers, retaining their COSM/COSN nomenclature. These identifiers remain the same 
# between different assemblies (GRCh37 and GRCh38). The legacy mutation identifier is 
# referred to as LEGACY_MUTATION_ID in the download files.

# Previously, each mutation at a specific genomic coordinate but on a different transcript
# had a unique COSM identifier. Now, all COSM identifiers at the same genomic location 
# have been collapsed into one representative COSM identifier. All previous COSM 
# identifiers are being maintained in order to enable tracking of existing mutations. 
# Mutations can be found in the website using legacy mutation identifiers, either by 
# entering the full identifier in the search bar, e.g. COSM476, or by entering the URL 
# directly, inserting the numerical part of the COSM, 
# e.g. https://cancer.sanger.ac.uk/cosmic/mutation/overview?id=476.

# TO DO: The GENOMIC_MUTATION_ID is not the preferred mutation identifier so don't do 
# filtering based on LEGACY_MUTATION_ID. Just check if there are still duplicated 
# LEGACY_MUTATION_ID after filtering based on GENOMIC_MUTATION_ID.

indG <- which(duplicated(ncv.df$GENOMIC_MUTATION_ID))
indL <- which(duplicated(ncv.df$LEGACY_MUTATION_ID))
x$GLsame <- ifelse(identical(indG, indL), TRUE, FALSE)
x$dupGid <- length(indG)
x$dupLid <- length(indL)
rm(indL)

# Mutations in a site is marked by a unique, sample-independent GENOMIC_MUTATION_ID.
# Check that entries with the same GENOMIC_MUTATION_ID do have the same position
ncv.df <- ncv.df[,c("GENOMIC_MUTATION_ID", "genome position")]
dupGid <- ncv.df$GENOMIC_MUTATION_ID[duplicated(ncv.df$GENOMIC_MUTATION_ID)]
dupGid.ind <- which(ncv.df$GENOMIC_MUTATION_ID%in%dupGid)
same.TF <- by(data=ncv.df$`genome position`[dupGid.ind],
              INDICES=ncv.df$GENOMIC_MUTATION_ID[dupGid.ind],
              FUN=function(x) length(unique(as.character(x)))==1
)

x$dupGidSamepos <- ifelse( all(same.TF), TRUE, FALSE )

# Drop entries with same GENOMIC_MUTATION_ID but diff position
GidDrop.ind <- dupGid.ind[!same.TF]
x$dupGidDiffpos <- length(GidDrop.ind)
x$DropSampidGid <- length(unique(which(SampDrop.TF), GidDrop.ind)) 

#rm( list=ls()[!ls()%in%c(x, out.dir, genome)] ); gc()

save(x, file=paste0(out.dir, "/CosmicNCV_GRCh", genome, "_explore.RData"))

# Extract sample dataset
ncv.df <- as.data.frame(ncv.df[1:sample.len,])
save(ncv.df, file=paste0(out.dir, "/CosmicNCV_GRCh", genome, "_", nrow(ncv.df), ".RData"))
rm(mut); gc()

#rm(list=ls()); gc()
