################################################################################
# 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
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
src.dir = paste0(wk.dir, "/z_ignore_git/out_mut_contact_Cp_plot")
out.dir = paste0(wk.dir, "/z_ignore_git/out_mut_contact_Cp_plot_sigsSummary_cor")

mutsig.file = paste0(data.dir, "/signal_mutSig/out_samplesForSignature/donorlist_signatureExposureRAW.csv")
### OTHER SETTINGS #############################################################
Cp.v = 1:21

mut.data.id = "donor_centric_PCAWG_Hg19"
sig.filter.id = "sigEperclimits_nosampfilter_ijmut"

mut.calcs = c("numWTSEQ", "Tmut", "Tmutnorm", "Nmsite", "Nmsitenorm", "TmutDIVNmsite")
mut.types = "All" #c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
mut.types = gsub(">", "To", mut.types, fixed=T)
mut.locs = "exon" #c("exon", "intron", "intergenic", "intron_intergenic", 
                  #  "exon_intron_intergenic_intergenic.excl")

# Get list of signatures
sig.df = read.csv(file=mutsig.file)[,-(1:7)]
mut.sigs = "RefSig.1" #c("RefSig.MMR1_RefSig.MMR2", colnames(sig.df))
rm(sig.df)

# Cp MEAN or MED (median)?
average.fnx = "MEAN"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
source(paste0(lib, "/doCorTest.R")) # Update deva copy
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
combi.df <- expand.grid(type=mut.types, loc=mut.locs, calc=mut.calcs, sig=mut.sigs, 
                        stringsAsFactors=F)
combi.len <- length(combi.df$type)

# Data for CSV

df <- foreach(c.ind=1:combi.len 
                   
  ) %do% {
  
  type <- combi.df$type[[c.ind]]
  loc <- combi.df$loc[[c.ind]]
  calc <- combi.df$calc[[c.ind]]
  sig <- combi.df$sig[[c.ind]]
  
  src.nme <- paste0(calc, "_", mut.data.id, "_", type, "_", sig, "_", loc, "_", sig.filter.id)
  
  # Load pairwisedifftest.RData
  load(paste0(src.dir, "/chrALL_", src.nme, "_pairwisedifftest.RData"))
  
  # Mean, median, count per Cp
  
  src.mx <- do.call("rbind", TEST$meanmed)
  Cp1To21.averages <- src.mx[as.character(Cp.v), average.fnx]
  
  # Correlation coefficient to represent trend of means
  
  COR <- doCorTest(xval=as.numeric(Cp.v), yval=Cp1To21.averages, 
                   alt="two.sided",method.names=c("pearson", "spearman", "kendall"), 
                   exactpval=F, out.dir=NULL, out.name=NULL)
  
  # Extract coefficients, p-values
  COR <- sapply(X=c("pear", "spea", "kend"), simplify=F, USE.NAMES=T, FUN=function(cor.meth){
    cor.vals <- COR[[cor.meth]]
    cor.vals <- c(coeff=unname(cor.vals$estimate), pval=cor.vals$p.value)
    return(cor.vals)
  })
  
  COR.and.N <- c(unlist(COR), list(min.N=min(src.mx[,"N"]), max.N=max(src.mx[,"N"])))
  
  return(COR.and.N)
  
}

df <- lapply(df, FUN=as.data.frame)
df <- do.call("rbind", df)
df <- cbind(combi.df, df)

df$pear.pval.BH <- p.adjust(p=df$pear.pval, method="BH")
df$spea.pval.BH <- p.adjust(p=df$spea.pval, method="BH")
df$kend.pval.BH <- p.adjust(p=df$kend.pval, method="BH")

save(df, file=paste0(out.dir, "/chrALL_", mut.data.id, "_", sig.filter.id, 
                "_Cp", average.fnx, "_cor_N.RData"))
write.csv(df, file=paste0(out.dir, "/chrALL_", mut.data.id, "_", sig.filter.id, 
                          "_Cp", average.fnx, "_cor_N.csv"), row.names=F)


# rm(list=ls()); gc()











