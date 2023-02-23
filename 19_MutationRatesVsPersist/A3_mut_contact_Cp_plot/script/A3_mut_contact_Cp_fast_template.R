################################################################################
# Per sig, calc, type, loc, generate df (with data from all chr), df.stat (with
# summary statistics of df including Cp=0 or all long-range contacts), calculate
# P-values comparing Cp distributions and calc vs. Cp, generate average calc
# values vs. Cp per signature (but mainly for nosampfilter). 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/19_MutationRatesVsPersist")
#src.dir = paste0(wk.dir, "/out_mut_contact_Cp_plotdata")
src.dir = paste0(wk.dir, "/out_mut_contact_Cp_plotdata_PCAWG_RefSig.1_nosampfilter")
out.dir = paste0(wk.dir, "/out_mut_contact_Cp_plot_fast")

mutsig.file = paste0(data.dir, "/signal_mutSig/out_samplesForSignature/donorlist_signatureExposureRAW.csv")
### OTHER SETTINGS #############################################################
Cp.v = 1:21

mut.data.id = "ijfnxmean_donor_centric_PCAWG_Hg19"
sig.filter.id = "sigEperclimits_nosampfilter_ijmut" #"sigEperclimits_1000rawInf_ijmut"

chrs = paste0("chr", c(1:22))
nCPU = 1
#mut.calcs = c("numWTSEQ", "Tmut", "Nmsite", "Tmutnorm", "Nmsitenorm", "TmutDIVNmsite")
mut.calcs = c("Tmutnorm", "Nmsitenorm", "TmutDIVNmsite")
mut.types = c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
mut.types = gsub(">", "To", mut.types, fixed=T)
mut.locs = c("exon", "intron", "intergenic", "intron_intergenic", "exon_intron_intergenic_intergenic.excl")

# Get list of signatures
sig.df = read.csv(file=mutsig.file)
mut.sigs = colnames(sig.df)[grepl("RefSig.", colnames(sig.df))]
mut.sigs = "RefSig.1" #c(mut.sigs, "RefSig.MMR1_RefSig.MMR2")
rm(sig.df)

combi.df.ind = INDREPLACE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))

library(Rmisc)
#source(paste0(lib, "/doVarTest_new.R")) # Update deva copy
source(paste0(lib, "/doCorTest.R")) # Update deva copy
source(paste0(lib, "/compareManyDist.R"))  # Update deva copy
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrs.len <- length(chrs)

combi.df <- expand.grid(loc=mut.locs, type=mut.types, calc=mut.calcs, sig=mut.sigs, 
                        stringsAsFactors=F)[combi.df.ind,]
combi.len <- length(combi.df$type)

foreach(c.ind=1:combi.len, .combine="rbind"
                   
) %do% {
  
  type <- combi.df$type[[c.ind]]
  loc <- combi.df$loc[[c.ind]]
  calc <- combi.df$calc[[c.ind]]
  sig <- combi.df$sig[[c.ind]]
  
  # Per calc, sig, type, loc, combine data from all chr -> df
  
  src.nme <- paste0(calc, "_", mut.data.id, "_", type, "_", sig, "_", loc, "_", sig.filter.id)
  
  toExport <- c("src.dir", "chrs", "src.nme", "calc", "sig", "type", "loc")
  #### PARALLEL EXECUTION #########
  df <- foreach(itr=isplitVector(1:chrs.len, chunks=nCPU), .inorder=F, .combine="rbind",
                .export=toExport, .noexport=ls()[!ls()%in%toExport]
                
  ) %op% {
    
    chunk <- sapply(itr, simplify=F, FUN=function(i){
      
      # Load IJ.MUT
      src.file <- paste0(src.dir, "/", chrs[[i]], "_", src.nme, ".RData")
      if( file.exists(src.file) ){
        load(src.file)
        return(IJ.MUT)
      } else {
        return(NULL)
      }
      
    })
    return( do.call("rbind", chunk))
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  if( is.null(df) ){
    message(paste0(src.nme, ": No data from any chr, skipped."))
    rm(df)
    next
  }
  
  df <- as.data.frame(df)
  df$Cp <- factor(as.character(df$Cp), levels=c(as.character(Cp.v)))
  
  # df summary statistics
  df.stat.Cp1To21 <- summarySE(data=df, measurevar="value", groupvars="Cp",
                               na.rm=T, .drop=F)
  
  # Add summary statistics for Cp=0 (all long-range contacts)
  df.stat.Cp0 <- summarySE(data=df, measurevar="value", groupvars=NULL, na.rm=T, .drop=F)[,-1]
  df.stat.Cp0 <- cbind.data.frame(Cp=factor("0", level="0"), df.stat.Cp0)
  
  C0To21.meds <- c(median(df$value, na.rm=T), 
                   aggregate(df$value, by=list(df$Cp), median, na.rm=T)$x)
  
  # Final df.stat
  df.stat <- cbind(calc=calc, sig=sig, type=type, loc=loc, 
                   rbind(df.stat.Cp0, df.stat.Cp1To21),
                   med=C0To21.meds)
  
  #save(df, file=paste0(out.dir, "/chrALL_", src.nme, ".RData"))
  save(df.stat, file=paste0(out.dir, "/chrALL_", src.nme, "_summ_stat.RData"))
  
  ## P-values
  
  # P-values, ANOVA/KW and correlation tests
  
  df <- na.omit(df)
  
  #try(doVarTest( xval=df$value, grp=df$Cp, out.dir=out.dir, out.name=paste0("chrALL_", src.nme) ))
  
  try(doCorTest( xval=as.numeric(as.character(df$Cp)), yval=df$value, alt="two.sided",
                 exactpval=F, out.dir=out.dir, out.name=paste0("chrALL_", src.nme) ))
  
  ## Add all values as Cp=0
  
  #df.Cp0 <- df
  #df.Cp0$Cp <- factor("0", levels="0")
  #df <- rbind(df.Cp0, df)
  
  try(compareManyDist( xval=df$value, grp=df$Cp, alt="two.sided", out.dir=out.dir, 
                       out.name=paste0("chrALL_", src.nme) ))
  rm(df)
  
} # combi.len foreach loop end


# rm(list=ls()); gc()