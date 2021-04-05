################################################################################
# Explore COSMIC non-coding mutation data with and without filtering.
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
genome = 38
CosmicNCVPath = paste0(data.dir, "/cosmic/GRCh", genome, "/CosmicNCV.GRCh", 
                       genome, ".tsv")
out.dir = paste0(wk.dir, "/out_explore")
### OTHER SETTINGS #############################################################
sample.len = 100000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
library(data.table)
source(paste0(wk.dir, "/lib/exploreData.R"))
source(paste0(wk.dir, "/lib/applyFilter1.R"))
source(paste0(wk.dir, "/lib/applyFilter2.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
X <- list()

ncv.df <- fread(file=CosmicNCVPath, header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)

#-------------------Initial exploration

X$orig <- exploreData(df=ncv.df)

#-------------------Filter1

# Take only mutations from WGS samples
ncv.df <- ncv.df[ncv.df$Whole_Exome=="n" & ncv.df$Whole_Genome_Reseq=="y",]

# Initial filtering based on sample and GENOMIC_MUTATION_ID
filter1 <- applyFilter1(ncv.df=ncv.df)
write.csv(x=filter1$SampleRefTable, file=paste0(out.dir, "/CosmicNCV_GRCh", genome, 
                                              "_SampleRefTable.csv"), row.names=FALSE)

ncv.df <- ncv.df[filter1$SAMP.TF & filter1$NOTREDUNDANT.TF,]
rm(filter1); gc()

#-------------------Filter2

filter2 <- applyFilter2(ncv.df=ncv.df)
X$dupGIDperSampAfterFilter1 <- sum(filter2$GMI.drop.TF )
X$dupPosperSampAfterFilter1 <- sum(filter2$POS.drop.TF)

ncv.df <- ncv.df[!(filter2$GMI.drop.TF | filter2$POS.drop.TF),]
rm(filter2); gc()

#-------------------Explore after Filter 1 and 2

X$Filtered1and2 <- exploreData(df=ncv.df)

#-------------------Extract sample dataset
set.seed(123)
ncv.df <- ncv.df[sample(x=1:nrow(ncv.df), size=sample.len),]
sample.len <- format(x=nrow(ncv.df), scientific=FALSE)
save(ncv.df, file=paste0(out.dir, "/CosmicNCV_GRCh", genome, "_", sample.len, ".RData"))

save(X, file=paste0(out.dir, "/CosmicNCV_GRCh", genome, "_explore.RData"))

# rm(list=ls()); gc()
