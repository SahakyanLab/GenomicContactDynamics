################################################################################
# Generate MUTBIN.DF containing mutation calculations for each mutation type,
# signature/s of interest, signature exposure threshold and mutation location.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Expands warnings
#options(warn=1)

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
src.dir = paste0(data.dir, "/cancer_somatic_mutation/out_filter")
BASECONT.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binBaseContent")
out.dir = paste0(wk.dir, "/out_mutCalcPerBin")

sigExposure.dir = paste0(data.dir, "/signal_mutSig/out_samplesForSignature")
chrLen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
data.id = "donor_centric_PCAWG" # "CosmicNCV", "donor_centric_ICGC" "donor_centric_PCAWG"
src.id = "Hg19" # "Hg19" | "hg38ToHg19"
bin.len = 40000L
# Metadata will contain total mutations and number of samples per MUTBIN.DF generated
metadata.id = "comprehensivefinal"

mut.v = "T>G" #c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

# Should match levels of ncv.df$location 
loc.v = c("exon", "intron", "intergenic", "intron_intergenic", 
          "exon_intron_intergenic_intergenic.excl")

# Directories in BASECONT.dir/ containing the BINKMER.MX
basecont.id.v = c(exon="maskfile1", intron="maskfile2", intergenic="maskfile3",
                  `intron_intergenic`="maskfile4", 
                  `exon_intron_intergenic_intergenic.excl`="maskfile0")
# BINKMER.MX affix, _<affix> except for ""
basecont.affix.v = c(exon="_hg19_Transcript_intron_hg19_intergenic_0updownTr", 
                     intron="_hg19_Transcript_exon_hg19_intergenic_0updownTr", 
                     intergenic="_hg19_Transcript_exon_hg19_Transcript_intron_hg19_intergenic_2000excluded",
                     `intron_intergenic`="_hg19_Transcript_exon_hg19_intergenic_2000excluded", 
                     `exon_intron_intergenic_intergenic.excl`="")

# Load signature exposure table to get signatures
sigExpRAW.df <- read.csv(file=paste0(sigExposure.dir, "/donorlist_signatureExposureRAW.csv"), 
                         header=T, stringsAsFactors=F) 
sigExpRAW.df <- sigExpRAW.df[,-c(1,3:7)]
sigExpPERC.df <- read.csv(file=paste0(sigExposure.dir, "/donorlist_signatureExposurePERCENT.csv"), 
                          header=T, stringsAsFactors=F)
sigExpPERC.df <- sigExpPERC.df[,-c(1,3:7)]
print(paste0("Signature exposure dfs has ", ncol(sigExpPERC.df[,-1]), " signatures..."), quote=FALSE)

# Define signatures
SIG.v = c("RefSig.MMR1_RefSig.MMR2", 
          colnames(sigExpPERC.df)[colnames(sigExpPERC.df)!="alt.ID"])

# Only consider samples/donor within this signature exposure percentage range (right closed interval).
sigEpLim.v = c("nosampfilter", "-1_0", "0_100", "5_100", "10_100", "30_100", "50_100", "1000rawInf") 
#"nosampfilter" - skip filtering based on signature exposure, use all samples in the dataset 

nCPU = 1L # Number of combinations of mut.v, SIG.v, sigEpLim.v and loc.v
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(wk.dir, "/lib/selectSamplesBasedOnSigExposure.R"))
source(paste0(wk.dir, "/lib/makeMUTBINDFperChrPerMUT.R"))
source(paste0(wk.dir, "/lib/generateMUTBINDF.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
mut.v <- unique(mut.v)
SIG.v <- unique(SIG.v)
sigEpLim.v <- unique(sigEpLim.v)
loc.v <- unique(loc.v)

if( any(!loc.v%in%names(basecont.id.v)) ){
  stop("Some locations not represented in basecont.id.v.")
}
if( any(!loc.v%in%names(basecont.affix.v)) ){
  stop("Some locations not represented in basecont.affix.v.")
}

chrlen.df <- read.delim(file=chrLen.file, header=T, stringsAsFactors=F)

for(mut in mut.v){
  
  # Generate combinations of mut, sig.v, sigEpLim.v, loc.v
  comb.mx <- data.matrix(expand.grid(1:length(mut), 1:length(SIG.v), 1:length(sigEpLim.v), 
                                     1:length(loc.v), stringsAsFactors=F))
  i.len <- nrow(comb.mx)
  
  toExport <- c("src.dir", "comb.mx", "mut", "SIG.v", "sigEpLim.v", 
                "loc.v", "basecont.id.v", "basecont.affix.v", "data.id", "src.id", "out.dir",
                "BASECONT.dir", "bin.len", "chrlen.df", "sigExpRAW.df", "sigExpPERC.df")
  
  #### PARALLEL EXECUTION #########
  metadata <- foreach(itr=isplitVector(1:i.len, chunks=nCPU), .inorder=T, 
                      .combine="rbind", .export=toExport, 
                      .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    
    # Load mutation data
    load(file=paste0(src.dir, "/", data.id, "_", src.id, "_final.RData"))
    #load(file=paste0(src.dir, "/", data.id, "_", src.id, "_final_100000.RData"))
    
    # Filter by mutation type
    if(mut=="All"){
      print("No filtering by mutation type.")
    } else if( mut%in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") ){
      ncv.df <- ncv.df[ncv.df$MUT==mut,]
    } else {
      stop(mut, ": Invalid mutation type notation...")
    }
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      
      v <- unname(comb.mx[i,])
      mut <- mut[ v[1] ]
      SIG <- SIG.v[ v[2] ]
      sigEpLim <- sigEpLim.v[ v[3] ]
      loc <- loc.v[ v[4] ]
      basecont.id <- basecont.id.v[loc]
      basecont.affix <- basecont.affix.v[loc]
      
      run.id <- paste0(data.id, "_", src.id, "_", mut, "_", SIG, "_", loc, 
                       "_sigEperclimits_", sigEpLim)
      print(paste0("RUN: ", run.id), quote=F)
      
      if( grepl(x=sigEpLim, pattern="_", fixed=T) | sigEpLim=="nosampfilter" ){
        sigExposure.df <- sigExpPERC.df
        print("Using percentage signature exposures...", quote=F)
      } else if( grepl(x=sigEpLim, pattern="raw", fixed=T) ){
        sigExposure.df <- sigExpRAW.df
        print("Using raw signature exposures...", quote=F)
      } else {
        stop(paste0("Invalid sigEpLim:", sigEpLim), quote=F)
      }
        
      return(
        generateMUTBINDF(ncv.df=ncv.df, data.id=data.id, src.id=src.id, out.dir=out.dir, 
                         basecont.dir=paste0(BASECONT.dir, "/", basecont.id), 
                         basecont.affix=basecont.affix, bin.len=bin.len, chrlen.df=chrlen.df, 
                         mut=mut, sigExposure.df=sigExposure.df, SIG=SIG, sigEpLim=sigEpLim, loc=loc)
      )
      
      #rm(v, mut, SIG, sigEpLim, loc, basecont.id, basecont.affix)
      
    })
    
    return( do.call("rbind", chunk) )
    #return( do.call("rbind.data.frame", chunk) )
    #return( do.call("rbind.data.frame", c(chunk, stringsAsFactors=F)) )
    #return(chunk)
    
  }
  
  ### END OF PARALLEL EXECUTION ###
  
  colnames(metadata) <- c("mut.id", "SIG.id", "loc.id", "sigEpLim.id", "Nmut", "Nsamp",
                          "meanNsamp", "medianNsamp", "sample")
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=T)
  write.csv(x=metadata, file=paste0(out.dir, "/", data.id, "_", src.id, "_", mut.id, "_", metadata.id, ".csv"),
            row.names=F)
  
} # mut.v for loop end

# rm(list=ls()); gc()