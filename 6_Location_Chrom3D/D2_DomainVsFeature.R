################################################################################
# Map features to each domain in DOMXYZR.DF. Add every feature to DOMXYZR.DF as 
# an extra column with values corresponding to their count per domain.
# deva, R/3.5.0-newgcc, gcc/4.9.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    repeat.dir = "/Users/ltamon/Database/ucsc_tables/hsa_RepeatMasker"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    repeat.dir  = "/t1-data/user/ltamon/Database/ucsc_tables/hsa_RepeatMasker"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    repeat.dir  = "/home/ltamon/Database/ucsc_tables/hsa_RepeatMasker"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "IMR90_LMNB1_GSE49341_hg19" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
domXYZR.dir = paste0(wk.dir, "/out_AddXYZR")
out.dir = paste0(wk.dir, "/out_DomainVsFeature")
### OTHER SETTINGS #############################################################
ploidy = "haploid"
out.name = paste0(model.id, "_", ploidy)
# Should correspond
featurefile.v = c(#paste0(repeat.dir , "/RepeatMasker_hg19")
                  #,
                  paste0(wk.dir, "/out_CpAsFeature/", out.name,
                         "_min2Mb_CpFeat.RData")
                  #,
                  #paste0(wk.dir, "/out_CpAsFeature/", out.name,
                  #       "_min05Mb_CpFeat.RData")
)
# Feature name; Should correspond with featurefile.v
feature.v = c(#"hg19Repeats"
              #,
              "min2Mb_CpFeat"
              #,
              #"min05Mb_CpFeat"
              )

featRobj = TRUE
# Separate files per chromosome?
featSepChr = FALSE
# If featSepChr = TRUE, directory should only contain the feature files with 
# "chr<N>" in filename
# If featSepChr = TRUE, enumerate chr files for feature 
chr.v = NULL
nCPU = 5L #~16G

feat.chr   = "chr"
feat.start = "start"
feat.end   = "end"
feat.name  = "Cp"

#feat.chr   = "genoName"
#feat.start = "genoStart"
#feat.end   = "genoEnd"
#feat.name  = "repName"

# Overlap parameter
type.olap="any"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
library(GenomicRanges)
source(paste0(lib, "/loadRData.R"))
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
### FUNCTION ###################################################################
DomainVsFeature <- function(
  out.dir = paste0(wk.dir, "/out_Chrom3D"),
  domXYZR.dir = paste0(wk.dir, "/out_Chrom3D"),
  out.name = "IMR90",
  feature = "hg19Repeats",
  featurefile = paste0(feat.dir, "/RepeatMasker_hg19"),
  featRobj = FALSE,
  # Separate files per chromosome?
  featSepChr = FALSE,
  # If featSepChr = TRUE, directory should only contain the feature files with 
  # "chr<N>" in filename
  # If featSepChr = TRUE, enumerate chr files for feature 
  chr.v = NULL,
  nCPU = 20L,
  feat.chr   = "genoName",
  feat.start = "genoStart",
  feat.end   = "genoEnd",
  feat.name  = "repName",
  # Overlap parameter
  type.olap="any" # within
){
  
  if(nCPU > 1){
    registerDoParallel(cores=nCPU)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", nCPU, " cores."), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  # DOMXYZR.DF
  DOMAIN.df <- loadRData(file=paste0(domXYZR.dir, "/", out.name, "_domXYZR.RData"))
  
  # Load feature file
  if(featRobj==FALSE & featSepChr==FALSE){
    FEATURE.df <- fread(file=featurefile, header=TRUE, data.table=FALSE,
                        stringsAsFactors=FALSE)
  } else if(featRobj==TRUE & featSepChr==FALSE){
    FEATURE.df <- loadRData(file=featurefile)
  }
  
  # Get common chromosomes of domainfile and featurefile
  if(featSepChr==FALSE){
    common.chr <- intersect( unique(DOMAIN.df$chr), FEATURE.df[[feat.chr]] )
  } else {
    common.chr <- intersect( unique(DOMAIN.df$chr), chr.v )
  }
  
  allfeature.uniq <- sort( unique(FEATURE.df[[feat.name]]) )
  
  common.chr.len <- length(common.chr)
  toExport <- c("common.chr", "featRobj", "featSepChr", "feat.dir", "feat.header",
                "DOMAIN.df", "FEATURE.df", "feat.chr","feat.start", "feat.end",
                "type.olap", "feat.name", "allfeature.uniq")
  
  #### PARALLEL EXECUTION #########
  
  DOMXYZRFEAT.DF <- foreach( itr=isplitVector(1:common.chr.len, chunks=nCPU),
                             .combine="rbind", .inorder=TRUE,
                             .export=toExport, .noexport=ls()[!ls()%in%toExport] 
  ) %op% {
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      chr <- common.chr[i]
      
      if(featRobj==FALSE & featSepChr==TRUE){
        FEATURE.df <- fread(file=paste0(feat.dir, "/", dir(feat.dir, pattern=chr)), 
                            header=feat.header, data.table=FALSE, stringsAsFactors=FALSE)
      } else if(featRobj==TRUE & featSepChr==TRUE){
        FEATURE.df <- loadRData(file=paste0(feat.dir, "/", dir(feat.dir, pattern=chr)))
      }
      
      domain.df <- DOMAIN.df[DOMAIN.df$chr==chr,]
      feature.df <- FEATURE.df[FEATURE.df[[feat.chr]]==chr,]
      
      # Query   <- domain
      # Subject <- feature
      olap.mx <- WhichOverlap(start.query=domain.df[["start"]], 
                              end.query=domain.df[["end"]], 
                              space.query=rep("a", nrow(domain.df)),
                              start.subject=feature.df[[feat.start]], 
                              end.subject=feature.df[[feat.end]], 
                              space.subject=rep("a",nrow(feature.df)),
                              maxgap=-1L, minoverlap=1L,
                              type=type.olap)
      
      # Per domain: vector of frequency per feature 
      by.domain <- by(data=olap.mx[,"subject"], #features
                      INDICES=domain.df[olap.mx[,"query"], "id"], #domains
                      FUN=function(f.ind){ #indices (feature.df)
                        feat <- feature.df[as.numeric(f.ind), feat.name]
                        #add NA as factor if present
                        feat <- addNA(x=feat, ifany=TRUE) 
                        feat.Freq <- as.data.frame( table(feat, useNA="ifany") )
                        feat.Freq[ match(x=allfeature.uniq, table=feat.Freq$feat),
                                   "Freq" ]
                      })
      
      # df with domains as rows and features as columns
      # Eeach Freq vector per domain from lst is a row
      df <- do.call("rbind", by.domain) 
      df <- data.frame(rownames(df), df, row.names=NULL,
                       stringsAsFactors=FALSE)
      colnames(df) <- c("id", allfeature.uniq)
      
      # all.x=TRUE to keep domains that don't have any overlapping repeats
      # These domains all have NA
      # merge() does not keep order of x
      merge(x=domain.df, y=df, by="id", all.x=TRUE)
      
    })
    return(do.call("rbind", chunk))
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  rm(allfeature.uniq); gc()
  
  rownames(DOMXYZRFEAT.DF) <- NULL
  save(DOMXYZRFEAT.DF, file=paste0(out.dir, "/", out.name, "_", feature, 
                                   "_domXYZRFeat.RData"))
  
}
################################################################################
DomainVsFeature <- cmpfun(DomainVsFeature, options=list(suppressUndefined=TRUE))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
feature.v.len <- length(feature.v)

if(feature.v.len!=length(featurefile.v)){
  stop("Checkpoint 1.")
}

for(i in 1:feature.v.len){
  
  DomainVsFeature(
    out.dir=out.dir,
    domXYZR.dir=domXYZR.dir,
    out.name=out.name,
    feature=feature.v[i],
    featurefile=featurefile.v[i],
    featRobj=featRobj,
    # Separate files per chromosome?
    featSepChr=featSepChr,
    # If featSepChr = TRUE, directory should only contain the feature files with 
    # "chr<N>" in filename
    # If featSepChr = TRUE, enumerate chr files for feature 
    chr.v=chr.v,
    nCPU=nCPU,
    feat.chr=feat.chr,
    feat.start=feat.start,
    feat.end=feat.end,
    feat.name=feat.name,
    # Overlap parameter
    type.olap=type.olap # within
  )
  
}

# rm(list=ls())



