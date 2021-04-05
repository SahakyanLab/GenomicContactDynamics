################################################################################
# Plot median and mean of mutation calculations per Cp. Calculate significant
# relative to Cp=1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/binmx/out_bindata_1perc_HiCNorm"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
mutCp.dir = paste0(wk.dir, "/out_mutCalcPerCp")
out.dir = paste0(wk.dir, "/out_mutCalcPerCp_sc")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
data.id = "CosmicNCV" # "CosmicNCV" | "donor_centric_PCAWG_as_is" | "donor_centric_PCAWG_sigEperc30_MMR12"
src.id = "hg38ToHg19" # "Hg19" | "hg38ToHg19"
mut.v = c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

calc.v = "Nmsitenorm" # c("Tmut", "Nmsite", "TmutDIVNmsite", "Nmsitenorm")
calc.id = "Nmsitenorm" # "Allcalc" | "Nmsitenorm"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/aggregateDF.R"))
source(paste0(wk.dir, "/lib/identifyAltHyp.R"))
source(paste0(wk.dir, "/lib/doMannWhitneyAndScatter.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(gcb, "_", data.id, "_", src.id, "_", calc.id)

p.lst <- NULL
for(mut in mut.v){
  
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=TRUE)
  
  mutcp.id <- paste0(gcb, "_", data.id, "_", src.id, "_", mut.id)
  load(file=paste0(mutCp.dir, "/", mutcp.id, "_mutCalcPerCp.RData"))
  rm(mutcp.id)
  MUTCP.DF <- MUTCP.DF[,c(calc.v, as.character(1:21))]
 
  #-------------------Boxplots
  x <- reshape2::melt(data=MUTCP.DF, id.vars=calc.v)
  rm(MUTCP.DF); gc()
  # Take only bins with Cp >= 1
  x <- x[x$value==1,colnames(x)!="value"]
  colnames(x)[colnames(x)=="variable"] <- "ind"
  x$ind <- as.numeric(x$ind)
  
  for(calc in calc.v){
    
    df <- aggregateDF(ind=x$ind, values=x[[calc]])
    id <- paste0(out.id, "_", mut.id, "_", calc)
    p.lst[[id]] <- doMannWhitneyAndScatter(x=x, df=df, calc=calc, plot.id=id)$p
    rm(df, id); gc()
    
  }  # calc.v for loop end
  
  print(paste0(mut, " done!"), quote=FALSE)
  rm(x); gc()
  
} # mut.v for loop end

mut.v.len <- length(mut.v)
calc.v.len <- length(calc.v)

if(calc.v.len>1){
  # For Allcalc
  p.arr <- ggarrange(plotlist=p.lst, nrow=mut.v.len, ncol=calc.v.len, legend=NULL)
  ggexport(p.arr, width=calc.v.len*10, height=mut.v.len*10,
           filename=paste0(out.dir, "/", out.id, "_scatplot.pdf"))
} else if(calc.v.len==1){
  # For single calc
  p.arr <- ggarrange(plotlist=p.lst, ncol=mut.v.len, nrow=calc.v.len, legend=NULL)
  ggexport(p.arr, height=calc.v.len*10, width=mut.v.len*10,
           filename=paste0(out.dir, "/", out.id, "_scatplot.pdf"))
} else {
  stop("Are you joking?!")
}

# rm(list=ls()); gc()



