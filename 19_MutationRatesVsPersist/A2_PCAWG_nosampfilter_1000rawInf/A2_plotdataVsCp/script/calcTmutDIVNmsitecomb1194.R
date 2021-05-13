################################################################################
# Plot mutation calculations vs. Cp (bin- and contact-wise) as boxplot and 
# scatterplot. Also, calculate significance of values across Cps relative to Cp=1.
# Generate plot per mutation type and signature. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Expands warnings
options(warn=1)

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/binmx/out_bindata_1perc_HiCNorm"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}

BASECONT.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binBaseContent")
mutbin.dir = paste0(wk.dir, "/out_mutCalcPerBin")
out.dir = paste0(wk.dir, "/out_plotdataVsCp")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_persistmx_ijonly")

chrLen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
# Contains all combinations of loc, mut, SIG, sigEpLim
comb.file = paste0(wk.dir, "/out_makecombfile/combfile_nosampfilter_1000rawInf")
### OTHER SETTINGS #############################################################
data.id = "donor_centric_PCAWG" # "CosmicNCV" | "donor_centric_PCAWG" | "donor_centric_ICGC"
src.id = "Hg19"                 # "Hg19" | "hg38ToHg19"
bin.len = 40000L
Cp.v = 1:21
# For contact-wise
gcb = "min2Mb" 
# Calculation in MUTBIN.DF
calc = "TmutDIVNmsite" 
# Row number of loc, mut, SIG, sigEpLim in comb.file
comb.num = 1194

# Central value to be used for scatterplot
# String of a function name e.g. "mean", "median"
aggregatefunx = "mean" 

nCPU = 1L # Contact-wise; chr
#-------------------
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
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
library(ggsci)
library(matrixStats)
library(reshape2)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/binToContact.R"))
source(paste0(wk.dir, "/lib/funxv.R"))
source(paste0(wk.dir, "/lib/getRelevantBins.R"))
source(paste0(wk.dir, "/lib/makeMUTBINDFtoBinMUTCPDFperMUT.R"))
source(paste0(wk.dir, "/lib/makeMUTBINDFtoContactMUTCPDFperMUT.R"))
source(paste0(wk.dir, "/lib/aggregateDF.R"))
source(paste0(wk.dir, "/lib/identifyAltHyp.R"))
source(paste0(wk.dir, "/lib/doMannWhitney.R"))
source(paste0(wk.dir, "/lib/makebp.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Get loc, mut, sigEpLim and SIG
comb <- readLines(con=comb.file)[comb.num]
comb <- strsplit(x=comb, split=";", fixed=T)[[1]]
if(length(comb)!=4){ sop("Invalid comb.") }
loc <- comb[4]
mut <- comb[2]
sigEpLim <- comb[3]
SIG <- comb[1]
rm(comb)

# Get total bins per chromosome
chrlen.df <- read.delim(file=chrLen.file, header=T, stringsAsFactors=F)
tot.bin.v <- ceiling(chrlen.df$length.bp/bin.len)
names(tot.bin.v) <- chrlen.df$chromosome
if( !identical(as.numeric(chrlen.df$bins.40kb), as.numeric(tot.bin.v)) ){
  stop("Checkpoint 1.")
}
rm(chrlen.df, chrLen.file)

# Function for combing values of two contacting bins
funx.v = c("ROWMEANS", "ROWSDS")
funx.id.v = c("mean", "sd")
wise.v <- c("b", paste("c", funx.id.v, sep=""))
wise.v.len <- length(wise.v)

mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=T)

#-------------------
# Load MUTBIN.DF
mutbin.id <- paste0(data.id, "_", src.id, "_", mut.id, "_", SIG, "_", loc, 
                    "_sigEperclimits_", sigEpLim)
mutbin.file <- paste0(mutbin.dir, "/", mutbin.id, "_mutCalcPerBin.RData")
if( file.exists(mutbin.file) ){
  load(file=mutbin.file)
} else {
  next
}

out.id <- paste0(gcb, "_", aggregatefunx, "_", calc, "_", mutbin.id)

#MUTBIN.DF <- MUTBIN.DF[MUTBIN.DF$chr%in%c("chr21", "chr22"),] # REMOVE

# Collect data for mut-sig plot
chr.str <- sort(unique(MUTBIN.DF$chr))
chr.str <- paste0(length(chr.str), "chrs_", paste(chr.str, collapse="."))

P.DF   <- list()
PFC.DF <- list()
pdf(file=paste0(out.dir, "/", out.id, "_bp.pdf"), height=10, width=10)
par(mfrow=c(2,2))

for(w in 1:wise.v.len){
  
  if(w==1){
    
    # Bin-wise
    x <- makeMUTBINDFtoBinMUTCPDFperMUT(binmx.dir=binmx.dir, 
                                        gcb=gcb, Cp.v=Cp.v, 
                                        bin.len=bin.len,
                                        MUTBIN.DF=MUTBIN.DF[,c("chr", "bin", calc)],
                                        basecont.dir=paste0(BASECONT.dir, "/", basecont.id.v[loc]), 
                                        basecont.affix=basecont.affix.v[loc])
  
    data.table::setnames(x=x, old=calc, new=wise.v[1], skip_absent=F)
    # Original number of Cp bins
    obin.len <- sum(!duplicated(x[,c("chr", "bins")]))
    # Take only bins with Cp>1
    x <- x[x$hasInd==1,]
    
  } else if(w==2){
    
    # Contact-wise
    rm(x, obin.len); gc()
    
    x <- makeMUTBINDFtoContactMUTCPDFperMUT(persist.dir=persist.dir,
                                            basecont.dir=paste0(BASECONT.dir, "/", basecont.id.v[loc]), 
                                            basecont.affix=basecont.affix.v[loc], 
                                            gcb=gcb, MUTBIN.DF=MUTBIN.DF[,c("chr", "bin", calc)],
                                            tot.bin.v=tot.bin.v, Cp.v=Cp.v,
                                            nCPU.chr=nCPU, funx.v=funx.v, funx.id.v=funx.id.v)
    
    obin.len <- nrow(x)
    data.table::setnames(x=x, old=funx.id.v[w-1], new=wise.v[w], skip_absent=F)
    rm(MUTBIN.DF); gc()
    
  } else {
    
    x$ind <- factor(as.character(x$ind), levels=as.character(Cp.v))
    data.table::setnames(x=x, old=funx.id.v[w-1], new=wise.v[w], skip_absent=F)
  
  } 
  
  val.id <- wise.v[w]
  
  # Take only bins with Cp>=1, bins overlapping with loc; this filtering 
  # should remove non-finite values.
  x <- x[!is.na(x$mutbin), colnames(x)!="hasInd"]
  if( any(!is.finite(x[[val.id]])) ){
    stop(paste0(mutbin.id, ": Non-finite ", val.id, " values." ))
  }
  
  # Count valid bins and contacts 
  vbin.len <- ifelse( wise.v[w]=="b", sum(!duplicated(x[,c("chr", "bins")])), nrow(x) )
  vbin.id <- paste0( wise.v[w], "_", vbin.len, "valid(", (vbin.len/obin.len)*100, 
                     "%)outof", obin.len, "\n", chr.str)
  x <- x[,!colnames(x)%in%c("chr","bins", "i", "j")]
  
  #-------------------Make boxplot
  if( !identical(levels(x$ind), as.character(Cp.v) ) ){
    stop(paste0(mutbin.id, ": Cp levels wrong." ))
  }
  
  # Make boxplot
  BINMUT <- makebp(df=x, calc=val.id, xlab="Cp", ylab=val.id, addjitter=F, 
                   plot.id=paste0(gcb, "_", mutbin.id, "_", calc, "\n", vbin.id))
  
  # x$ind left as factor for aggregateDF()
  DF.id <- paste0(mutbin.id, ".", val.id)
  P.DF[[DF.id]] <- aggregateDF(ind=x$ind, values=x[[val.id]], FUN=aggregatefunx)
  P.DF[[DF.id]] <- cbind.data.frame(mut.id=rep(mut.id), SIG.id=rep(SIG), loc.id=rep(loc), 
                                    sigEpLim.id=rep(sigEpLim), P.DF[[DF.id]], 
                                    percbinmut=as.numeric(BINMUT$percbinmut[ P.DF[[DF.id]]$ind ]),
                                    binPerCp=as.numeric(BINMUT$binPerCp[ P.DF[[DF.id]]$ind ]),
                                    wise=wise.v[w], stringsAsFactors=F)
  
  rm(BINMUT, vbin.len, vbin.id)
  
  # Mann-Whitney
  x$ind <- as.character(x$ind)
  P.DF[[DF.id]] <- doMannWhitney(x=x, df=P.DF[[DF.id]], calc=val.id)
  rownames(P.DF[[DF.id]]) <- P.DF[[DF.id]]$ind
  
  # Fold change version relative to Cp=1
  PFC.DF[[DF.id]] <- P.DF[[DF.id]]
  PFC.DF[[DF.id]]$values <- log2(P.DF[[DF.id]]$values/P.DF[[DF.id]]["1","values"])
  
  print(paste0(wise.v[w], " ", sigEpLim, " ", SIG, " ", loc, " done!"), quote=F)
  x[[val.id]] <- NULL
  rm(val.id); gc()
  
} # wise.v.len for loop end

dev.off()

P.DF <- do.call("rbind.data.frame", c(P.DF, stringsAsFactors=F))
PFC.DF <- do.call("rbind.data.frame", c(PFC.DF, stringsAsFactors=F))
rownames(P.DF) <- rownames(PFC.DF) <- NULL

save(P.DF, file=paste0(out.dir, "/", out.id, "_scatplot.RData"))
save(PFC.DF, file=paste0(out.dir, "/", out.id, "_FC_scatplot.RData"))

# rm(list=ls()); gc()