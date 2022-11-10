################################################################################
# Assign mean/median Cp to Phylo-HMRF hg19 converted regions
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) 
options(warn=1)

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/27_EvolutionaryAnalyses")
src.dir = paste0(data.dir, "/Phylo-HMRF/out_combine_splitPerChr") # 0-based
ijLO.dir = paste0(wk.dir, "/out_liftover")
out.dir = paste0(wk.dir, "/out_assignCp")
### OTHER SETTINGS #############################################################
src.id = "genome_state_Phylo-HMRF_mapping_contact50K_norm"
src.header = T
chrs = "chrCHRREPLACE"
LOchain = "hg38ToHg19"
LOwidth.min.bp = 30000 # Take only regions maintaining length (50000) after liftover
# as well as output regions >= LOwidth.min.bp (latter is relevant for regions
# broken down into multiple regions after conversion i.e. 1 long plus few short ones)
# Ideally, choose value > half of original length to get 1 output region for 1
# input region.
coord.system = 1 
consensus.FUN = "mean"
# PERSIST.MX
gcb = "min2Mb"
bin.len = 40000 
nCPU = 2 # Contacts
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/getConsensusCpOfPairs.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chrs){
  
  # 
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  
  ij.mx <- cbind(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
  rm(PERSIST.MX)
  
  #
  load(paste0(ijLO.dir, "/", chr, "_", LOchain, "_LOwidth.min.bp", LOwidth.min.bp, 
              "_", src.id, ".RData"))
  ijLO.mx <- data.matrix(ijLO.df[,-1])
  rm(ijLO.df)

  #
  
  consensusCp <- getConsensusCpOfPairs(bin.len=bin.len, nCPU=nCPU, 
                                       ij.mx=ij.mx, # [i, j, Cp]  
                                       coord.system=coord.system,
                                       p.coord.mx=ijLO.mx, # 1-based CS, [start.a, end.a, start.b, end.b]
                                       consensus.FUN=consensus.FUN)
  
  ij.pair.len <- length(ijLO.mx[,1])
  if( !identical(as.integer(consensusCp[,"pair.ind"]), 1:ij.pair.len) ){
    rm(consensusCp)
    stop(paste(chr, ": consensusCp order not consistent with ijLO.df/mx."))
  }

  consensusCp <- consensusCp[,"value", drop=F]
  save(consensusCp, file=paste0(out.dir, "/", chr, "_", LOchain, "_LOwidth.min.bp", 
                                LOwidth.min.bp, "_", src.id, "_", consensus.FUN, "Cp.RData"))
  
  print(paste0(chr, ": done!"), quote=F)
  
  rm(ijLO.mx, consensusCp, chr)
  gc()
  
}

# rm(list=ls()); gc()