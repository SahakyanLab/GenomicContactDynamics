################################################################################
### FUNCTION ###################################################################
applyFilter1 <- function(ncv.df=ncv.df){
  
  warning("Make sure ncv.df was filtered for only Whole_Genome_Reseq or Whole_Exome 
          mutations.")
  
  x <- list()
  
  #-------------------Get unique pairs of ID_SAMPLE-Sample name
  # There can be multiple ID_SAMPLE for a sample if the same sample has been entered into 
  # the database multiple times from different papers or if it is unclear whether the same 
  # sample name between different publications is indeed the same sample.
  
  # Remove ambiguity. Take only unique pairs of ID_SAMPLE-Sample name
  custom.id <- paste0(ncv.df$ID_tumour, ncv.df$`Primary site`, ncv.df$`Site subtype 1`, 
                      ncv.df$`Site subtype 2`, ncv.df$`Site subtype 3`,
                      ncv.df$`Primary histology`, ncv.df$`Histology subtype 1`, 
                      ncv.df$`Histology subtype 2`, ncv.df$`Histology subtype 3`)
  tmp1 <- by(data=paste0(ncv.df$ID_SAMPLE, custom.id),
             INDICES=ncv.df$`Sample name`, FUN=function(x){
               length(unique(x))
             })
  tmp1 <- tmp1[tmp1==1]
  tmp2 <- by(data=paste0(ncv.df$`Sample name`, custom.id),
             INDICES=ncv.df$ID_SAMPLE, FUN=function(x){
               length(unique(x))
             })
  tmp2 <- tmp2[tmp2==1]
  rm(custom.id)
  
  x$SAMP.TF <- ncv.df$`Sample name`%in%names(tmp1) & ncv.df$ID_SAMPLE%in%names(tmp2)
  rm(tmp1, tmp2)
  tmp <- ncv.df[x$SAMP.TF,c("Sample name", "ID_SAMPLE")]
  tmp <- tmp[ !duplicated(tmp), ]
  x$SampleRefTable <- tmp
  rm(tmp); gc()
  
  #-------------------
  # Remove redundant entries due to mutations in a sample being uploaded multiple 
  # times from different studies.
  v <- c("Sample name", "ID_SAMPLE", "GENOMIC_MUTATION_ID", "zygosity", "genome position",
         "Mutation somatic status", "WT_SEQ", "MUT_SEQ", "SNP", "HGVSG")
  x$NOTREDUNDANT.TF <- !duplicated(ncv.df[,v]) 
 
  return(x)

}

# rm(list=ls()); gc()

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

# TO DO: The GENOMIC_MUTATION_ID is the preferred mutation identifier so don't do 
# filtering based on LEGACY_MUTATION_ID. Just check if there are still duplicated 
# LEGACY_MUTATION_ID after filtering based on GENOMIC_MUTATION_ID.