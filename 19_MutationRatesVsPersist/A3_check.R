################################################################################
# Make sure that per sample, each site has only 1 or no mutation. 
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
src.dir = paste0(wk.dir, "/out_filter")
out.dir = paste0(wk.dir, "/out_check")
### OTHER SETTINGS #############################################################
src.id = "hg38ToHg19" # "hg38ToHg19" | "Hg19"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
#load(file=paste0(src.dir, "/CosmicNCV_", src.id, "_final_2000.RData"))
load(file=paste0(src.dir, "/CosmicNCV_", src.id, "_final.RData"))

# ID_SAMPLE and Sample name
# There can be multiple ids, if the same sample has been entered into the 
# database multiple times from different papers. -- Check for this cases

# Also, there can be multiple sample ids linked to the same sample name if it 
# is unclear whether the same sample name between different publications is 
# indeed the same sample.
custom.id <- paste(ncv.df$`Primary site`, ncv.df$`Site subtype 1`, 
                   ncv.df$`Site subtype 2`, ncv.df$`Site subtype 3`,
                   ncv.df$`Primary histology`, ncv.df$`Histology subtype 1`,
                   ncv.df$`Histology subtype 2`, ncv.df$`Histology subtype 3`)

# Trim down dataset
ncv.df <- ncv.df[, c("Sample name", "ID_SAMPLE", "GENOMIC_MUTATION_ID", "genome position")]

x <- list()
# Make sure 1 mut per pos per samp
x$onemutperpospersamp <- c(pos_ID_SAMPLE=NA, pos_Sample_name=NA, pos_custom.id=NA)
tmp <- paste(custom.id, ncv.df$`genome position`)
x$onemutperpospersamp["pos_custom.id"] <- ifelse( any(duplicated(tmp)), FALSE, TRUE )
rm(tmp, custom.id); gc()
tmp <- paste(ncv.df$ID_SAMPLE, ncv.df$`genome position`)
x$onemutperpospersamp["pos_ID_SAMPLE"] <- ifelse( any(duplicated(tmp)), FALSE, TRUE )
rm(tmp)
tmp <- paste(ncv.df$`Sample name`, ncv.df$`genome position`)
x$onemutperpospersamp["pos_Sample_name"] <- ifelse( any(duplicated(tmp)), FALSE, TRUE )
rm(tmp); gc()

# Entries of the same mutation type may have different GENOMIC_MUTATION_ID if from
# different samples. Check that entries with the same GENOMIC_MUTATION_ID do have the 
# same position
ncv.df <- ncv.df[,c("GENOMIC_MUTATION_ID", "genome position")]
dupGid <- ncv.df$GENOMIC_MUTATION_ID[duplicated(ncv.df$GENOMIC_MUTATION_ID)]
indG <- which(ncv.df$GENOMIC_MUTATION_ID%in%dupGid)
same.TF <- by(data=ncv.df$`genome position`[indG],
              INDICES=ncv.df$GENOMIC_MUTATION_ID[indG],
              FUN=function(x) length(unique(x))==1
              )
x$dupGidSamepos <- ifelse( all(same.TF), TRUE, FALSE )
rm(same.TF, dupGid, indG, ncv.df); gc()

# Given that there is 1 mut per pos per samp based on checks above and that 
# entries with the same GENOMIC_MUTATION_ID have the same position, all entries
# with the same GENOMIC_MUTATION_ID should be from different samples.

save(x, file=paste0(out.dir, "/CosmicNCV_", src.id, "_check.RData"))

# rm(list=ls()); gc()
