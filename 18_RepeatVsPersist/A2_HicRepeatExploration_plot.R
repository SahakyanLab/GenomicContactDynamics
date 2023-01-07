################################################################################
#  Make density of MINREP.MX values. This was used to generate the distribution
# of sumrep values.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

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

wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
rep.group = "subfamALL"
metric = "sumrep"
minelm.dir = out.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group, "_", metric)
### OTHER SETTINGS #############################################################
suffix = "min2Mb"
chr.v = paste0("chr", c(1:22, "X")) #"chr21" #"chrCHRREPLACE" 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chr.v){

  load(file=paste0(minelm.dir, "/", chr, "_MinElm_", suffix, ".RData"))
  
  pdf(file=paste0(out.dir, "/", chr, "_MinElm_", suffix, "_plot.pdf"),
      height=10, width=10)
  
  plot(density(as.numeric(log10(MINELM.MX[,-1])), na.rm=T), 
       main=paste0(rep.group, "_", metric, "_valueslog10transformed"))
  
  dev.off()
  
  print(paste0(chr, " done!"), quote=F)
  
  rm(MINELM.MX)
  gc()
  
}

# rm(list=ls()); gc()